classdef (Abstract) CropOp < handle
    methods (Abstract)
        run(obj, crop_dir, varargin)
    end
end
