classdef AffineRegistration < Registration
    % AffineRegistration  Frame-by-frame affine registration.
    %
    %   Supports translation, rotation, scaling, and shearing. Uses
    %   MATLAB's imregtform / imwarp (Image Processing Toolbox) with a
    %   per-frame affine2d transform estimated against the reference image.
    %
    %   Usage:
    %       reg = AffineRegistration(movie_in, reference_img);
    %       reg = reg.run();
    %
    %   Properties:
    %       imreg_type      - 'monomodal' (default) or 'multimodal'
    %                         passed to imregconfig()
    %       optimizer_overrides - struct of optimizer property overrides,
    %                             applied on top of imregconfig defaults.
    %                             E.g. struct('MaximumIterations', 200)

    properties
        % Re-instantiate inherited properties
        data_processed Movie

        imreg_type          char    = 'monomodal'
        optimizer_overrides struct  = struct()
    end

    methods (Access = protected)
        function [obj, movie_result, operation] = registration_alg(obj, data_raw)
            arguments
                obj   AffineRegistration
                data_raw Movie
            end

            % initialize
            movie_result = data_raw;
            movie = double(data_raw.stack);
            ref_img = double(obj.reference_img);
            if isempty(ref_img)
                ref_img = mean(movie, 3, 'omitmissing');
            end
            obj.reference_img = ref_img;

            ref_ref = imref2d(size(ref_img));
            n_frames = size(movie, 3);

            try
                tforms = obj.operation.tforms;
            catch
                tforms = [];
            end

            if isempty(tforms)
                disp('Computing affine transforms..')

                [optimizer, metric] = imregconfig(obj.imreg_type);

                % Apply any user-specified optimizer overrides
                fields = fieldnames(obj.optimizer_overrides);
                for i = 1:numel(fields)
                    optimizer.(fields{i}) = obj.optimizer_overrides.(fields{i});
                end

                tforms = cell(n_frames, 1);
                for f = 1:n_frames
                    frame = movie(:, :, f);
                    tforms{f} = imregtform(frame, ref_img, 'affine', optimizer, metric);
                    movie(:, :, f) = imwarp(frame, tforms{f}, 'OutputView', ref_ref);
                end

            else
                disp('Applying input affine transforms..')

                for f = 1:n_frames
                    frame = movie(:, :, f);
                    movie(:, :, f) = imwarp(frame, tforms{f}, 'OutputView', ref_ref);
                end
            end

            movie_result.stack = movie;
            operation.tforms = tforms;
        end
    end

    methods
        % Constructor
        function obj = AffineRegistration(movie_in, reference_img)
            arguments
                movie_in        = ''
                reference_img double = []
            end
            obj = obj@Registration(movie_in, reference_img);

            obj.init.method = 'Affine Registration';

            obj.data_processed = obj.data_raw;
            obj.operation      = [];
        end
    end
end
