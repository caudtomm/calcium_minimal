classdef RigidRegistration < Registration
    properties
        % Re-instantiate inherited properties
        data_processed Movie
    end

    methods (Access = protected)
        function [obj,movie_result,operation] = registration_alg(obj, data_raw)
            arguments
                obj RigidRegistration
                data_raw Movie
            end

            % initialize
            movie_result = data_raw;
            movie = data_raw.stack;
            ref_img = obj.reference_img;
            if isempty(ref_img); ref_img = mean(movie,3,'omitmissing'); end
            obj.reference_img = ref_img;
            try
                shifts = obj.operation.shifts;
            catch
                shifts = [];
            end
            
            if isempty(shifts)
                disp('Computing shifts..')
                [shift_x, shift_y] = register_frames(movie, ref_img);
%                 shift_x = mean(shift_x) * ones(size(shift_x));
%                 shift_y = mean(shift_y) * ones(size(shift_y));
            else
                disp('Applying input shifts..')
                shift_x = shifts(:,1);
                shift_y = shifts(:,2);
            end
            movie = shift_data(movie,shift_x,shift_y);
            
            movie_result.stack = movie;
            operation.shifts = [squeeze(shift_x),squeeze(shift_y)];
        end
    end
    
    methods
        % Constructor
        function obj = RigidRegistration(movie_in, reference_img)
            arguments
                movie_in = ''
                reference_img double = []
            end
            obj = obj@Registration(movie_in,reference_img);

            % definition of init structure
            obj.init.method = 'Rigid Registration';

            % initialization
            obj.data_processed = obj.data_raw;
            obj.operation = [];
        end

        
    end
end