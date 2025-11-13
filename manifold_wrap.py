from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat
import mat73, h5py, itertools, pandas as pd
from tqdm import tqdm
import time

def _try_mat73(p):
    try:
        return mat73.loadmat(p)
    except Exception:
        return None

def _try_loadmat(p):
    try:
        return loadmat(p, squeeze_me=True, struct_as_record=False)
    except Exception:
        return None

def _extract_cell(d, key=None):
    if key and key in d:
        v = d[key]
    else:
        ks = [k for k in d.keys() if not k.startswith("__")]
        if len(ks) != 1:
            raise ValueError("ambiguous keys: " + str(ks))
        v = d[ks[0]]
    if isinstance(v, (list, tuple)):
        return list(v)
    if isinstance(v, np.ndarray) and v.dtype == object:
        return [v.item(i) if v.ndim == 1 else v[i] for i in range(v.size)]
    if isinstance(v, np.ndarray) and v.ndim >= 2 and v.dtype != object:
        return [v]
    raise TypeError("not a cell-like object")

def normalize_subj(subj):
    if isinstance(subj, str):
        return subj
    if isinstance(subj, np.ndarray):
        if subj.dtype.kind in ("U", "S"):
            return str(subj.squeeze())
        # char stored as uint16?
        if subj.dtype.kind in ("u", "i"):
            arr = subj.squeeze().tolist()
            try:
                return "".join(chr(c) for c in arr)
            except Exception:
                return str(arr)
        return str(subj.squeeze())
    return str(subj)

def normalize_stims(stims):
    if stims is None:
        return []
    if isinstance(stims, str):
        return [stims]
    if isinstance(stims, (list, tuple)):
        return [str(s) for s in stims]
    if isinstance(stims, np.ndarray):
        if stims.dtype == object:
            out = []
            for s in stims.ravel():
                if isinstance(s, str):
                    out.append(s)
                elif isinstance(s, np.ndarray):
                    out.append(str(s.squeeze()))
                else:
                    out.append(str(s))
            return out
        else:
            return [str(x) for x in stims.ravel().tolist()]
    return [str(stims)]

def load_file_struct(p, cell_key="data"):
    d = _try_loadmat(p)
    if d is None:
        d = _try_mat73(p)
        if d is None:
            with h5py.File(p, "r"):
                raise RuntimeError("HDF5 v7.3 not parsed; provide key")

    # data
    if cell_key and cell_key in d:
        data = d[cell_key]
    else:
        ks = [k for k in d.keys() if not k.startswith("__")]
        if len(ks) != 1:
            raise ValueError("ambiguous keys: " + str(ks))
        data = d[ks[0]]

    # turn data into list of arrays (like before)
    if isinstance(data, (list, tuple)):
        data_list = list(data)
    elif isinstance(data, np.ndarray) and data.dtype == object:
        data_list = [data.item(i) if data.ndim == 1 else data[i] for i in range(data.size)]
    elif isinstance(data, np.ndarray) and data.ndim >= 2 and data.dtype != object:
        data_list = [data]
    else:
        raise TypeError("not a cell-like object")

    subj = normalize_subj(d.get("subjID", ""))
    stims = normalize_stims(d.get("stims", []))

    return data_list, subj, stims

def save_result(out_dir, base, r1, r2, outm_list,
                src_i, idx_i, src_j, idx_j,
                subj_i, stims_i, subj_j, stims_j,
                save_mode="mat", csv_float_fmt="%.9g"):
    out_dir = Path(out_dir)
    pre_cell = np.empty((1, len(r1)), dtype=object)
    for t, arr in enumerate(r1):
        pre_cell[0, t] = np.asarray(arr)

    # normalize stims for savemat (cellstr)
    stims_i_arr = np.array([str(s) for s in stims_i], dtype=object).reshape(1, -1)
    stims_j_arr = np.array([str(s) for s in stims_j], dtype=object).reshape(1, -1)

    if save_mode == "csv":
        if hasattr(r2, "to_csv"):
            r2.to_csv(out_dir / f"{base}_metrics.csv",
                      index=True,
                      float_format=csv_float_fmt)
        savemat(out_dir / f"{base}.mat",
                {
                    "preprocessed": pre_cell,
                    "alignment": outm_list,
                    "src_i": np.array(src_i),
                    "idx_i": np.array(idx_i),
                    "src_j": np.array(src_j),
                    "idx_j": np.array(idx_j),
                    "subj_i": np.array(subj_i),
                    "subj_j": np.array(subj_j),
                    "stims_i": stims_i_arr,
                    "stims_j": stims_j_arr,
                },
                do_compression=True)
        return

    if hasattr(r2, "to_numpy"):
        met_cols = np.array(list(map(str, r2.columns)), dtype=object).reshape(1, -1)
        if r2.index.dtype == "object":
            met_idx = np.array(list(map(str, r2.index.tolist())), dtype=object).reshape(-1, 1)
        else:
            met_idx = np.asarray(r2.index.to_numpy()).reshape(-1, 1)
        met_vals = r2.to_numpy()
        m_struct = {
            "columns": met_cols,
            "index": met_idx,
            "values": met_vals,
            "index_name": np.array([str(r2.index.name) if r2.index.name is not None else ""], dtype=object),
        }
    else:
        m_struct = {
            "columns": np.array([], dtype=object).reshape(1, 0),
            "index": np.array([], dtype=float).reshape(0, 1),
            "values": np.asarray(r2),
        }

    savemat(out_dir / f"{base}.mat",
            {
                "preprocessed": pre_cell,
                "metrics": m_struct,
                "alignment": outm_list,
                "src_i": np.array(src_i),
                "idx_i": np.array(idx_i),
                "src_j": np.array(src_j),
                "idx_j": np.array(idx_j),
                "subj_i": np.array(subj_i),
                "subj_j": np.array(subj_j),
                "stims_i": stims_i_arr,
                "stims_j": stims_j_arr,
            },
            do_compression=True)

def process_input_path(in_path,
                       out_dir=None,
                       n_repetitions=150,
                       n_points=70,
                       cell_key="data",
                       save_mode="mat",
                       csv_float_fmt="%.9g"):
    import sys
    print(sys.executable)
    print(sys.version)

    from glue import preprocess as glue_pre
    from glue.contrib import glue_analysis_dataframe

    def preprocess(manifolds, idx):
        return glue_pre.downsample_manifolds(
            manifolds,
            n_points=n_points,
            seed=idx,
        )

    def analyze(manifolds, idx):
        return glue_analysis_dataframe(
            manifolds=manifolds,
            indices=(idx,),
            indices_name=["rep_num"],
            analysis_type="FIRST_VERSUS_REST",
            return_matrix=True,
            shuffle=True,
            gaussianize=False,
            bias=True,
            seed=idx,
        )

    in_path = Path(in_path)

    if in_path.is_dir():
        mat_files = sorted(in_path.glob("*.mat"))
        if out_dir is None:
            out_dir = in_path / "results"
    else:
        mat_files = [in_path]
        if out_dir is None:
            out_dir = in_path.parent / "results"

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mf = []
    src = []
    meta_subj = []
    meta_stims = []
    for p in mat_files:
        data_list, subj, stims = load_file_struct(p.as_posix(), cell_key)
        for idx_arr, arr in enumerate(data_list):
            a = np.asarray(arr)
            if a.ndim != 2 or a.dtype.kind != "f":
                a = a.astype(np.float64)
            mf.append(a)
            src.append((p.name, idx_arr))
            meta_subj.append(subj)
            meta_stims.append(stims)

    print("manifolds:", len(mf))
    pairs = list(itertools.combinations(range(len(mf)), 2))
    for i_f in mf:
        print("mf sz:", i_f.shape, flush=True)        

    for k, (ii, jj) in enumerate(tqdm(pairs, desc="pairs")):
        mfs = [mf[ii], mf[jj]]
        if any(m.shape[0] == 0 for m in mfs) or any(m.shape[1] < n_points for m in mfs):
            continue
        res = []
        outm_list = []
        for i_rep in range(n_repetitions):
            r1 = preprocess(mfs, i_rep)
            r2, outm = analyze(r1, i_rep)
            res.append(r2)
            outm_list.append(outm)
        r2_all = pd.concat(res, ignore_index=False)
        base = f"p_{ii:05d}_{jj:05d}"
        save_result(
            out_dir,
            base,
            r1,
            r2_all,
            outm_list,
            src[ii][0],
            src[ii][1],
            src[jj][0],
            src[jj][1],
            meta_subj[ii],
            meta_stims[ii],
            meta_subj[jj],
            meta_stims[jj],
            save_mode=save_mode,
            csv_float_fmt=csv_float_fmt,
        )


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="input .mat file OR directory with .mat files")
    p.add_argument("--outdir", required=False, help="output directory (default: <input>/results)")
    p.add_argument("--reps", type=int, default=50)
    p.add_argument("--n_points", type=int, default=70)
    p.add_argument("--save_mode", choices=["mat", "csv"], default="mat")
    args = p.parse_args()

    start = time.time()
    process_input_path(
        in_path=args.input,
        out_dir=args.outdir,
        n_repetitions=args.reps,
        n_points=args.n_points,
        save_mode=args.save_mode,
    )
    print("Time elapsed:", time.time() - start)
