classdef Snippet < Movie
    properties
        parentpath
        frameinterval double
    end

    methods (Access=private)
        function obj = setPath(obj)
            if isempty(obj.path)
                warning(['No parent reference found. Please make sure you ' ...
                    'know where this Snippet is from!'])
            end

            obj.parentpath = obj.path;
            obj.path = [];  % optionally store the output path for this 
                            % Snippet here, later.
        end
    end

    methods
        function obj = Snippet(src,frameinterval)
            % arg 'src' is a filename, numerical matrix or Movie object
            arguments
                src
                frameinterval double = [1]
            end            
            
            srcinput = src;
            if isa(src,'Movie')
                srcinput = src.stack;
            end
            
            % superclass constructor
            obj = obj@Movie(srcinput);

            % cut stack
            obj.stack = obj.stack(:,:,frameinterval);
            obj = obj.checkAgain(); % update frameavg and other properties
            obj.frameinterval = frameinterval(:);
            
            % set other poperties
            if isa(src,'Movie')
                obj.path = src.path;    % this will be cut&pasted into 
                                        % parentpath with method setPath()
                obj.badperiods = src.badperiods;
                if ~isempty(src.scanimage_meta)
                    obj.scanimage_meta = src.scanimage_meta(obj.frameinterval);
                end
            else
                if ~isempty(obj.scanimage_meta)
                    % in case it was loaded from a mat file containing a Movie object
                    obj.scanimage_meta = obj.scanimage_meta(obj.frameinterval);
                end
            end
            obj = setPath(obj); % set parent and output paths
            % update bad periods to new interval (make sure to update 
            % property 'frameinterval' before calling this)
            obj.badperiods = obj.update_badperiods();
        end

        function newbp = update_badperiods(obj)
            % this (1) retains only the badperiods pertaining to the 
            % Snippet frames, (2) ensures ordinal numbering corresponds to
            % the Snippet stack dim 3 coordinates
            
            bp = obj.badperiods;
            spacer = 1; % first col is assumed to be a homogeneous tag 
                        % (same value throughout, ex. trial num)

            if isempty(bp); newbp = bp; return; end
            
            % expand periods
            bplogical = convertPeriods(bp(:,spacer+[1:2]),true);
            
            % turn frame interval to logical
            frinterv = false(max(obj.frameinterval),1);
            frinterv(obj.frameinterval) = true;

            % equalize array lengths
            if length(bplogical) > length(frinterv)
                extlen = length(bplogical) - length(frinterv);
                ext = false(extlen,1);
                frinterv = [frinterv; ext];
            elseif length(bplogical) < length(frinterv)
                extlen = length(frinterv) - length(bplogical);
                ext = false(extlen,1);
                bplogical = [bplogical; ext];
            end

            % retain relevant badperiods only and convert to period matrix
            newbp = convertPeriods(bplogical(frinterv));
            if isempty(newbp)
                newbp = [];
            elseif spacer>0
                nbp = size(newbp,1);
                firstcol = repelem(bp(1,1),nbp);
                newbp = [firstcol(:), newbp];
            end
        end

        % save method (##store the output path to obj.path)
        function FileOut = save(obj, newfname, newpath, type) % ## ----- arguments should be the same for all Movie subclasses
            arguments
                obj
                newfname char
                newpath char = pwd
                type char = 'mat'
            end
            type = char(type);

            % define new fname (no ext!)
            % obj.path is assumed to be empty!
            obj.path = getFileNameSpecs('');
            if ~isempty(newfname)
                obj.path.fname = newfname;
            else
                basename = '';
                if ~isempty(obj.parentpath) && ~isempty(obj.parentpath.fname)
                    basename = [obj.parentpath.fname,'_'];
                end
                basename = [basename, '_Snippet#'];

                i = 1;
                while true
                    obj.path.fname = [basename,num2str(i)];
                    files = dir(fullfiletol(newpath,[obj.path.fname,'.*']));
                    if isempty(files)
                        break
                    else
                        i = i+1;
                    end
                end
            end

            % run inherited save function
            FileOut = save@Movie(obj, newpath, type);
        end
    end
end