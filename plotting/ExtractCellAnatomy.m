classdef ExtractCellAnatomy
% ExtractCellAnatomy  Load and crop per-trial image stacks for a single cell.
%
%   Loads fish1.mat, retrieves anatomy_imgs and localcorr_imgs
%   [H x W x T_total], crops both to a neighbourhood around the cell ROI,
%   and returns the binary ROI mask in crop coordinates.
%
%   Usage:
%       eca = ExtractCellAnatomy(traces, cell_idx);
%       eca = ExtractCellAnatomy(traces, cell_idx, 'neighborhood', 30);
%       [anat, lcorr, roi_mask, roi_label] = eca.getCroppedStack();
%
%   Inputs:
%       traces     ActivityTraces  — carries subject_locations, ROImap,
%                                    goodNeuron_IDs
%       cell_idx   scalar double   — 1-based index into goodNeuron_IDs
%
%   Options:
%       neighborhood  (default 20)  pixel padding around ROI bounding box
%
%   Outputs:
%       anat_imgs  [H_c x W_c x T_total]  cropped anatomy stack (all trials)
%       lcorr_imgs [H_c x W_c x T_total]  cropped local-corr stack, or []
%       roi_mask   [H_c x W_c]            logical ROI mask in crop coords
%       roi_label  scalar                 ROI label used in traces.ROImap
%
%   All outputs are [] / NaN if fish1.mat is missing or anatomy_imgs is empty.
%   lcorr_imgs is [] if localcorr_imgs is empty in the Subject file.
%
%   Note: populate anatomy_imgs by running subject.retrieve_trial_anatomies().

    properties
        traces      ActivityTraces
        cell_idx    double
        neighborhood double = 20
    end

    methods

        function obj = ExtractCellAnatomy(traces, cell_idx, options)
            arguments
                traces      ActivityTraces
                cell_idx    (1,1) double {mustBePositive, mustBeInteger}
                options.neighborhood (1,1) double {mustBePositive} = 20
            end
            obj.traces       = traces;
            obj.cell_idx     = cell_idx;
            obj.neighborhood = options.neighborhood;
        end

        function snip = loadTrialMovieSnippet(obj, trial_num, pre_interval)
        % loadTrialMovieSnippet  Load a full trial movie, compute dF/F, return a Snippet.
        %
        %   snip = loadTrialMovieSnippet(trial_num, pre_interval)
        %
        %   trial_num     scalar double  — 1-based trial number
        %   pre_interval  double array   — frame indices passed to Snippet
        %
        %   Returns [] if the movie file is not found or traces_src is unset.

            snip = [];

            loc   = obj.traces.subject_locations;
            fpath = fullfiletol(loc.subject_datapath, loc.traces_src);
            fname = [loc.subject_ID, '_', num2str(trial_num, '%05d'), '_*'];
            files = dir(fullfiletol(fpath, fname));

            if isempty(files)
                warning('ExtractCellAnatomy:movieNotFound', ...
                    'No movie file found for trial %d:\n  %s', trial_num, ...
                    fullfiletol(fpath, fname));
                return
            end

            raw_movie = robust_io('load', ...
                fullfiletol(files(1).folder, files(1).name), 'movie').movie;
            proc = BasicMovieProcessor('dff', raw_movie);
            evalc('proc = proc.run();');   % suppress per-call console noise
            dff = proc.data_processed;

            % Clamp pre_interval to valid frame range
            valid_fr = pre_interval(pre_interval >= 1 & pre_interval <= dff.nfr);
            if isempty(valid_fr)
                warning('ExtractCellAnatomy:intervalOutOfRange', ...
                    'Frame interval [%d %d] out of movie range [1 %d] for trial %d.', ...
                    pre_interval(1), pre_interval(end), dff.nfr, trial_num);
                return
            end
            snip = Snippet(dff, valid_fr);
        end

        function [r1, r2, c1, c2, roi_mask, roi_label] = getCropBounds(obj)
        % getCropBounds  Bounding box around the ROI in full-frame coordinates.
        %
        %   [r1, r2, c1, c2, roi_mask, roi_label] = getCropBounds()
        %
        %   Uses traces.ROImap directly — no fish1.mat loading needed.

            if obj.cell_idx > numel(obj.traces.goodNeuron_IDs)
                error('ExtractCellAnatomy:cellIdxOutOfRange', ...
                    'cell_idx=%d exceeds goodNeuron_IDs length (%d).', ...
                    obj.cell_idx, numel(obj.traces.goodNeuron_IDs));
            end
            roi_label     = obj.traces.goodNeuron_IDs(obj.cell_idx);
            roi_mask_full = (obj.traces.ROImap == roi_label);
            if ~any(roi_mask_full(:))
                error('ExtractCellAnatomy:roiNotFound', ...
                    'ROI label %d (cell_idx=%d) not found in ROImap.', ...
                    roi_label, obj.cell_idx);
            end
            [H, W] = size(obj.traces.ROImap);
            [r1, r2, c1, c2, roi_mask] = obj.computeCropBounds(roi_mask_full, H, W);
        end

        function [anat_imgs, lcorr_imgs, roi_mask, roi_label] = getCroppedStack(obj)
        % getCroppedStack  Load Subject, crop anatomy and localcorr stacks.
        %
        %   [anat_imgs, lcorr_imgs, roi_mask, roi_label] = getCroppedStack()

            anat_imgs = []; lcorr_imgs = []; roi_mask = []; roi_label = NaN;

            %% -- Load Subject ------------------------------------------------
            fish_path = fullfiletol( ...
                obj.traces.subject_locations.subject_datapath, 'fish1.mat');

            if ~exist(fish_path, 'file')
                warning('ExtractCellAnatomy:fileNotFound', ...
                    'Subject file not found — skipping:\n  %s', fish_path);
                return
            end

            subject = robust_io('load', fish_path).fish1;

            %% -- Validate anatomy images -------------------------------------
            anat = subject.anatomy_imgs;
            if isempty(anat)
                warning('ExtractCellAnatomy:noAnatomyImgs', ...
                    ['anatomy_imgs is empty — skipping.\n', ...
                     'Run subject.retrieve_trial_anatomies() and re-save.\n', ...
                     'File: %s'], fish_path);
                return
            end
            [H, W, ~] = size(anat);

            %% -- Locate ROI and compute crop bounds --------------------------
            if obj.cell_idx > numel(obj.traces.goodNeuron_IDs)
                error('ExtractCellAnatomy:cellIdxOutOfRange', ...
                    'cell_idx=%d exceeds goodNeuron_IDs length (%d).', ...
                    obj.cell_idx, numel(obj.traces.goodNeuron_IDs));
            end
            roi_label     = obj.traces.goodNeuron_IDs(obj.cell_idx);
            roi_mask_full = (obj.traces.ROImap == roi_label);

            if ~any(roi_mask_full(:))
                error('ExtractCellAnatomy:roiNotFound', ...
                    'ROI label %d (cell_idx=%d) not found in ROImap.', ...
                    roi_label, obj.cell_idx);
            end

            [r1, r2, c1, c2, roi_mask] = obj.computeCropBounds(roi_mask_full, H, W);

            %% -- Crop image stacks ------------------------------------------
            anat_imgs = cropStack(anat, r1, r2, c1, c2);

            lcorr = subject.localcorr_imgs;
            if ~isempty(lcorr)
                lcorr_imgs = cropStack(lcorr, r1, r2, c1, c2);
            end
        end

    end % public methods

    % ------------------------------------------------------------------
    methods (Access = private)

        function [r1, r2, c1, c2, roi_mask] = computeCropBounds(obj, roi_mask_full, H, W)
        % Bounding box of roi_mask_full expanded by obj.neighborhood, clamped
        % to image dims.  Also returns roi_mask in the cropped coordinate frame.
            row_coords = find(any(roi_mask_full, 2));
            col_coords = find(any(roi_mask_full, 1));
            nb = obj.neighborhood;
            r1 = max(1, min(row_coords) - nb);
            r2 = min(H, max(row_coords) + nb);
            c1 = max(1, min(col_coords) - nb);
            c2 = min(W, max(col_coords) + nb);
            roi_mask = roi_mask_full(r1:r2, c1:c2);
        end

    end % private methods

end % classdef


% ======================================================================
% File-private helper
% ======================================================================

function stack_out = cropStack(stack, r1, r2, c1, c2)
% Apply spatial crop [r1:r2, c1:c2] to a [H x W] or [H x W x T] array.
    if ndims(stack) == 2 %#ok<ISMAT>
        stack_out = stack(r1:r2, c1:c2);
    else
        stack_out = stack(r1:r2, c1:c2, :);
    end
end
