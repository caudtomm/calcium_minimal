"""
run_slicetca.py  —  Minimal sliceTCA workflow for neural activity tensors.

Reads per-subject *_input.mat files written by SliceTCA_Analysis.saveInputFiles
(MATLAB v7.3 / HDF5 format, tensor shape [trials x neurons x time]), then:

  1. Per-neuron min-max normalization (no temporal filtering)
  2. Grid search  — finds the number of components per slice type that
                    minimises reconstruction loss
  3. decompose    — full-resolution fit at the optimal ranks
  4. invariance   — post-hoc orthogonalisation of each slice type

Results are written as *_slicetca.mat files (scipy v5 format, loadable by
MATLAB's built-in load()).

Usage
-----
    python run_slicetca.py  INPUT_DIR  OUTPUT_DIR  [options]

    python run_slicetca.py  data/slicetca_inputs  data/slicetca_results \\
        --max_rank 6  --sample_size 3  --max_iter 15000  --seed 42

Component convention
--------------------
For input tensor shape (T=trials, N=neurons, S=time_bins), sliceTCA with
number_components=[r0, r1, r2] produces a 3-element list of partitions,
each partition being a 2-element list [scores, weights]:

    components[0][0]  shape (r0, T)     — trial-slice   scores  (one value per trial)
    components[0][1]  shape (r0, N, S)  — trial-slice   weights (neurons×time pattern)
    components[1][0]  shape (r1, N)     — neuron-slice  scores  (one value per neuron)
    components[1][1]  shape (r1, T, S)  — neuron-slice  weights (trials×time pattern)
    components[2][0]  shape (r2, S)     — time-slice    scores  (one value per time bin)
    components[2][1]  shape (r2, T, N)  — time-slice    weights (trials×neurons pattern)

Saved as components_{k}_scores and components_{k}_weights in each output .mat file.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
import traceback
from pathlib import Path
from typing import List

import numpy as np
import scipy.io
import torch

import mat73
import slicetca
from slicetca import decompose, grid_search, invariance


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _load_mat(path: str) -> dict:
    """Load a MATLAB .mat file regardless of version (v5 or v7.3/HDF5)."""
    try:
        return mat73.loadmat(path)
    except Exception:
        pass
    import scipy.io as sio
    return sio.loadmat(path, squeeze_me=True, struct_as_record=False)


def _str_from_mat(value) -> str:
    """Convert a MATLAB-loaded string to a plain Python str."""
    if isinstance(value, str):
        return value
    if isinstance(value, np.ndarray):
        v = value.squeeze()
        if v.dtype.kind in ("U", "S"):
            return str(v)
        # char arrays stored as uint16
        if v.dtype.kind in ("u", "i") and v.ndim == 1:
            try:
                return "".join(chr(c) for c in v.tolist())
            except Exception:
                return str(v.tolist())
    return str(value)


def _labels_from_mat(value) -> List[str]:
    """Convert a MATLAB cell-of-strings to a Python list of str."""
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple)):
        return [_str_from_mat(v) for v in value]
    if isinstance(value, np.ndarray):
        if value.dtype == object:
            return [_str_from_mat(v) for v in value.ravel()]
        if value.dtype.kind in ("U", "S"):
            return [str(s) for s in value.ravel()]
    return []


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

def normalize_per_neuron(data: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """Per-neuron min-max normalization across all trials and time points.

    Parameters
    ----------
    data : np.ndarray, shape (trials, neurons, time)

    Returns
    -------
    np.ndarray, same shape, each neuron rescaled to [0, 1].
    """
    d_min = data.min(axis=(0, 2), keepdims=True)   # (1, neurons, 1)
    d_max = data.max(axis=(0, 2), keepdims=True)   # (1, neurons, 1)
    return (data - d_min) / (d_max - d_min + eps)


# ---------------------------------------------------------------------------
# Grid search helpers
# ---------------------------------------------------------------------------

def run_grid_search(
    data_tensor: torch.Tensor,
    max_rank: int,
    min_rank: int,
    sample_size: int,
    processes_grid: int,
    processes_sample: int,
    seed: int,
    decompose_kwargs: dict,
) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    """Run sliceTCA grid search over all rank combinations.

    Returns
    -------
    loss_grid  : ndarray, shape (max_rank-min_rank,) * 3 + (sample_size,)
    seed_grid  : ndarray, same shape
    min_ranks  : list of ints, the lower bound used (same for all dimensions)
    """
    n_dims    = data_tensor.ndim                        # 3 for our tensors
    max_ranks = [max_rank] * n_dims
    min_ranks = [min_rank] * n_dims

    loss_grid, seed_grid = grid_search(
        data_tensor,
        max_ranks=max_ranks,
        min_ranks=min_ranks,
        sample_size=sample_size,
        processes_grid=processes_grid,
        processes_sample=processes_sample,
        seed=seed,
        **decompose_kwargs,
    )
    return loss_grid, seed_grid, min_ranks


def best_ranks_from_grid(
    loss_grid: np.ndarray,
    min_ranks: List[int],
) -> List[int]:
    """Return the rank combination with the lowest mean loss across seeds.

    Parameters
    ----------
    loss_grid : (..., sample_size) ndarray from grid_search
    min_ranks : offset added to each grid index
    """
    mean_loss = loss_grid.mean(axis=-1)                 # average over seeds
    best_idx  = np.unravel_index(np.argmin(mean_loss), mean_loss.shape)
    return [int(idx + offset) for idx, offset in zip(best_idx, min_ranks)]


# ---------------------------------------------------------------------------
# Per-subject pipeline
# ---------------------------------------------------------------------------

def process_subject(mat_path: str, args: argparse.Namespace) -> str:
    """Full sliceTCA pipeline for a single subject .mat file.

    Returns the path to the saved result file.
    """
    print(f"\n{'='*64}")
    print(f"Processing: {Path(mat_path).name}")

    # -- Load ----------------------------------------------------------------
    mat        = _load_mat(mat_path)
    data_raw   = np.array(mat["data"],      dtype=np.float64)  # (trials, N, T)
    trial_labs = _labels_from_mat(mat.get("trial_labels"))
    t_ax       = np.array(mat["t"],         dtype=np.float64).ravel()
    subject_id = _str_from_mat(mat.get("subject_id", ""))
    group      = _str_from_mat(mat.get("group",      ""))
    framerate  = float(np.squeeze(mat.get("framerate", 1.0)))

    T, N, S = data_raw.shape
    print(f"  Subject : {subject_id}   Group : {group}")
    print(f"  Tensor  : {T} trials × {N} neurons × {S} time bins")

    # Guard against NaNs/Infs before normalisation
    n_bad = (~np.isfinite(data_raw)).sum()
    if n_bad > 0:
        print(f"  WARNING : {n_bad} non-finite values → replaced with 0")
        np.nan_to_num(data_raw, nan=0.0, posinf=0.0, neginf=0.0, copy=False)

    # -- 1. Normalisation ----------------------------------------------------
    data_norm   = normalize_per_neuron(data_raw).astype(np.float32)
    device      = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  Device  : {device}")
    data_tensor = torch.from_numpy(data_norm).to(device)

    # Shared kwargs forwarded to every decompose() call inside grid_search
    decompose_kwargs = dict(
        positive      = args.positive,
        learning_rate = args.learning_rate,
        max_iter      = args.max_iter_grid,
        min_std       = args.min_std,
    )

    # -- 2. Grid search ------------------------------------------------------
    print(f"\n  [grid_search]  max_rank={args.max_rank}  "
          f"min_rank={args.min_rank}  sample_size={args.sample_size}")

    loss_grid, seed_grid, min_ranks = run_grid_search(
        data_tensor,
        max_rank        = args.max_rank,
        min_rank        = args.min_rank,
        sample_size     = args.sample_size,
        processes_grid  = args.processes_grid,
        processes_sample= args.processes_sample,
        seed            = args.seed,
        decompose_kwargs= decompose_kwargs,
    )

    best_ranks = best_ranks_from_grid(loss_grid, min_ranks)
    print(f"  Best ranks : {best_ranks}  "
          f"(trial-slice, neuron-slice, time-slice)")

    # Safeguard: if grid selects all zeros, fall back to [1, 1, 1]
    if all(r == 0 for r in best_ranks):
        print("  WARNING : all-zero ranks selected — defaulting to [1, 1, 1]")
        best_ranks = [1, 1, 1]

    # -- 3. Final decomposition ----------------------------------------------
    print(f"\n  [decompose]  ranks={best_ranks}  max_iter={args.max_iter}")
    components, model = decompose(
        data_norm,
        number_components = best_ranks,
        positive          = args.positive,
        learning_rate     = args.learning_rate,
        max_iter          = args.max_iter,
        min_std           = args.min_std,
        seed              = args.seed,
        progress_bar      = True,
    )

    # -- 4. Invariance optimisation ------------------------------------------
    print("  [invariance]")
    invariance(model)
    components = model.get_components(numpy=True)

    # -- Save ----------------------------------------------------------------
    stem     = Path(mat_path).stem.replace("_input", "")
    out_path = os.path.join(args.output_dir, f"{stem}_slicetca.mat")
    pt_path  = os.path.join(args.output_dir, f"{stem}_model.pt")

    # PyTorch model — saved separately (full object, reloadable with torch.load)
    torch.save(model, pt_path)
    print(f"  Model  → {pt_path}")

    # Reconstruction tensor (normalised scale, same shape as input)
    try:
        with torch.no_grad():
            reconstruction = model().cpu().numpy().astype(np.float32)
    except Exception as exc:
        print(f"  WARNING: could not compute reconstruction ({exc})")
        reconstruction = np.zeros(0, dtype=np.float32)

    results: dict = {
        "subject_id"    : subject_id,
        "group"         : group,
        "best_ranks"    : np.array(best_ranks, dtype=np.int32),
        "loss_grid"     : loss_grid.astype(np.float32),
        "losses"        : np.array(model.losses, dtype=np.float32),
        "reconstruction": reconstruction,   # (trials, neurons, time) normalised
        "t"             : t_ax,
        "framerate"     : framerate,
        "trial_labels"  : np.array(trial_labs) if trial_labs else np.zeros(0),
    }

    # Each partition k: [scores, weights]
    #   k=0 trial-slice:  scores (r0,T)    weights (r0,N,S)
    #   k=1 neuron-slice: scores (r1,N)    weights (r1,T,S)
    #   k=2 time-slice:   scores (r2,S)    weights (r2,T,N)
    for k, comps_k in enumerate(components):
        if comps_k and len(comps_k) >= 2:
            results[f"components_{k}_scores"]  = np.asarray(comps_k[0], dtype=np.float32)
            results[f"components_{k}_weights"] = np.asarray(comps_k[1], dtype=np.float32)
        else:
            results[f"components_{k}_scores"]  = np.zeros(0, dtype=np.float32)
            results[f"components_{k}_weights"] = np.zeros(0, dtype=np.float32)

    scipy.io.savemat(out_path, results)
    print(f"  Result → {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run sliceTCA on per-subject neural tensors.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input_dir",  help="Directory with *_input.mat files, or a single *_input.mat file")
    p.add_argument("output_dir", help="Directory for *_slicetca.mat results")

    gs = p.add_argument_group("Grid search")
    gs.add_argument("--max_rank",        type=int,   default=5,
                    help="Maximum rank per slice type")
    gs.add_argument("--min_rank",        type=int,   default=0,
                    help="Minimum rank per slice type")
    gs.add_argument("--sample_size",     type=int,   default=3,
                    help="Random seeds per grid point")
    gs.add_argument("--max_iter_grid",   type=int,   default=3000,
                    help="Max iterations per grid-search fit")
    gs.add_argument("--processes_grid",  type=int,   default=1,
                    help="Parallel processes across rank combinations")
    gs.add_argument("--processes_sample",type=int,   default=1,
                    help="Parallel processes across seeds")

    dc = p.add_argument_group("Decomposition")
    dc.add_argument("--max_iter",        type=int,   default=10000,
                    help="Max iterations for the final decompose() call")
    dc.add_argument("--learning_rate",   type=float, default=0.005,
                    help="Adam learning rate")
    dc.add_argument("--min_std",         type=float, default=1e-3,
                    help="Convergence threshold (loss std over iter_std window)")
    dc.add_argument("--positive",        action="store_true",
                    help="Constrain all components to non-negative values")
    dc.add_argument("--seed",            type=int,   default=7,
                    help="Base random seed")

    return p


def main() -> None:
    args = build_parser().parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # 'fork' (Linux default) + CUDA = deadlock; force 'spawn' when parallelising
    if args.processes_grid > 1 or args.processes_sample > 1:
        torch.multiprocessing.set_start_method("spawn", force=True)

    if os.path.isfile(args.input_dir):
        input_files = [args.input_dir]
    else:
        input_files = sorted(glob.glob(os.path.join(args.input_dir, "*_input.mat")))
    if not input_files:
        print(f"ERROR: no *_input.mat files in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(input_files)} subject file(s).")

    succeeded, failed = [], []
    for mat_path in input_files:
        try:
            out = process_subject(mat_path, args)
            succeeded.append(out)
        except Exception as exc:
            print(f"\nERROR processing {Path(mat_path).name}: {exc}", file=sys.stderr)
            traceback.print_exc()
            failed.append(mat_path)

    print(f"\n{'='*64}")
    print(f"Done.  {len(succeeded)}/{len(input_files)} subjects succeeded.")
    if failed:
        print("Failed:")
        for p in failed:
            print(f"  {p}")


if __name__ == "__main__":
    main()
