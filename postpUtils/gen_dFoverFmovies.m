classdef gen_dFoverFmovies < MovieProcessing
% gen_dFoverFmovies  Pre-compute per-trial dF/F movies and save to disk.
%
% Subclass of MovieProcessing that computes dF/F for a single raw trial
% movie following the same bad-period handling used in Registration:
% NaN/bad frames are stripped before F0 estimation so they cannot skew the
% baseline, then reinserted as NaN frames in the output.
%
% Single-trial usage (MovieProcessing API):
%   proc = gen_dFoverFmovies(raw_movie);
%   proc = proc.run();
%   dff_movie = proc.data_processed;
%
% Batch usage (also the compiled-binary entry point):
%   gen_dFoverFmovies.runForDirectory(movie_dir, outpath)
%
%   movie_dir  full path to the directory containing raw trial .mat files
%   outpath    directory where output .mat files will be written
%
% Each output file is named identically to its input trial file and
% contains a single variable 'movie' holding the full-frame dF/F Movie
% object.  Existing files are silently overwritten.

    properties
        data_processed  Movie  = Movie('')
        operation              = struct()
    end

    % ------------------------------------------------------------------
    methods

        function obj = gen_dFoverFmovies(movie_in)
            arguments
                movie_in = ''
            end
            obj = obj@MovieProcessing(movie_in);
            obj.init.method = 'gen_dFoverFmovies';
        end

        % --------------------------------------------------------------
        function obj = run(obj)
        % run  Compute dF/F for obj.data_raw; result in obj.data_processed.
            evalc('obj.disp_runheader()');

            movie = obj.data_raw;

            % Strip bad-period frames so they do not bias F0 estimation.
            rbp   = BasicMovieProcessor('remove_badperiods', movie);
            evalc('rbp = rbp.run();');
            clean = rbp.data_processed;

            % Compute per-pixel dF/F on clean frames only.
            dff_proc = BasicMovieProcessor('dff', clean);
            evalc('dff_proc = dff_proc.run();');
            dff_clean          = dff_proc.data_processed;
            obj.operation      = dff_proc.operation;
            obj.data_processed = dff_clean;
            evalc('obj = obj.tailsequence(''dff'');');

            % Reinsert NaN frames at original bad-period positions.
            abp = BasicMovieProcessor('add_badperiods_as_nans', obj.data_processed);
            evalc('abp = abp.run();');
            obj.data_processed = abp.data_processed;
        end

    end % methods

    % ------------------------------------------------------------------
    methods (Static)

        function runForDirectory(movie_dir, outpath)
        % runForDirectory  Batch entry point: process all trial movies in a directory.
        %
        %   gen_dFoverFmovies.runForDirectory(movie_dir, outpath)
        %
        %   Iterates over all .mat files in movie_dir in sorted order,
        %   computes full-frame dF/F movies, and saves each result to
        %   outpath under the same filename as the input.

            if ~exist(outpath, 'dir')
                mkdir(outpath);
            end

            files = dir(fullfile(movie_dir, 'TC_*.mat'));
            files = files(~[files.isdir]);
            [~, order] = sort({files.name});
            files = files(order);
            nfiles = numel(files);

            fprintf('Source  : %s  |  %d files\n', movie_dir, nfiles);
            fprintf('Output  : %s\n', outpath);

            for k = 1 : nfiles
                in_file = fullfile(files(k).folder, files(k).name);
                fprintf('\n[%d/%d] %s\n', k, nfiles, files(k).name);

                raw_movie = robust_io('load', in_file, 'movie').movie;

                proc = gen_dFoverFmovies(raw_movie);
                proc = proc.run();

                [~, base] = fileparts(files(k).name);
                proc.data_processed.save(outpath, 'mat', base, true);
            end

            fprintf('\nDone. %d files processed.\n', nfiles);
        end

    end % static methods

end % classdef
