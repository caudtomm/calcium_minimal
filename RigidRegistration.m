classdef RigidRegistration < Registration
    properties
    end
    
    methods
        % Constructor
        function obj = RigidRegistration()
            obj.method = 'Rigid';
        end

        function registration_alg()
            % alignment across planes ---- TEMP
            %if ~exist('ref_img', 'var') || isempty(ref_img); ref_img = mean(movie(:,:,40:300),3); end
            if ~exist('ref_img', 'var') || isempty(ref_img); ref_img = movie(:,:,1); end
            %movie = alginWithinTrial(movie, ref_img);
            if ~exist('shifts','var') || isempty(shifts)
                [shift_x, shift_y] = register_frames(movie, ref_img);
                shift_x = mean(shift_x) * ones(size(shift_x));
                shift_y = mean(shift_y) * ones(size(shift_y));
            else
                shift_x = shifts(:,1);
                shift_y = shifts(:,2);
            end
            movie = shift_data(movie,shift_x,shift_y);
        end
    end
end