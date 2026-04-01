classdef StackH5
    % StackH5
    % A lightweight, disk-backed proxy for a 3-D image stack stored in HDF5.
    % Goal: behave like a double HxWxF array for common usage patterns
    % (A(:,:,k), size(A), mean(A,3), mean(A,"all"), std(A,[],"all")),
    % while keeping RAM usage tiny and allowing huge movies.
    %
    % Notes
    % - Under the hood, data are kept as uint8 in HDF5 (lossless).
    % - All reads return double to match existing code that expects doubles.
    % - Only contiguous 1:K style indexing is supported along each dim.
    %   (e.g., A(:,:,1:10) OK; A(:,:,[1 3 7]) not supported -> stream with a loop.)
    % - fliplr/flipud are logical flags; the file is not rewritten.

    properties (SetAccess=private)
        f char              % HDF5 file path
        d char = '/stack'   % dataset path inside HDF5
        sz (1,3) double     % [H W F] size of the stack
        flipx logical = false % optional logical left-right flip flag
        flipy logical = false % optional logical up-down flip flag
    end

    methods
        function obj = StackH5(h5file,dset,sz,flipx,flipy)
            % Constructor
            % h5file : path to HDF5 file
            % dset   : dataset name in HDF5 (default '/stack')
            % sz     : [H W F]
            % flipx  : optional initial horizontal flip flag
            % flipy  : optional initial vertical flip flag
            obj.f = h5file;
            if nargin>=2 && ~isempty(dset), obj.d = dset; end
            obj.sz = sz;
            if nargin>=4, obj.flipx = logical(flipx); end
            if nargin>=5, obj.flipy = logical(flipy); end
        end

        function varargout = size(obj, dim)
            % size(A) / size(A,dim) / [m,n,p,...] = size(A)
            sz = obj.sz;

            if nargin > 1
                % size(A, dim)
                if dim <= numel(sz)
                    s = sz(dim);
                else
                    s = 1;
                end
                varargout = {s};
                return
            end

            if nargout <= 1
                % size(A) -> [H W F]
                varargout = {sz};
            else
                % [m,n,p,...] = size(A)
                % pad with 1s if caller asks for more dims than we have
                out = [sz, ones(1, nargout - numel(sz))];
                varargout = num2cell(out(1:nargout));
            end
        end

        function a = subsref(obj,S)
            % Provide A(...), A.prop and nested refs support.
            switch S(1).type
                case '.'
                    % Pass-through for property/method access.
                    a = builtin('subsref',obj,S);
                case '()'
                    % Handle array-style indexing on disk.
                    idx = S(1).subs;
                    if isscalar(idx) && isequal(idx{1},':')
                        % A(:) or A(:) followed by more subsref.
                        a = obj.read([1 1 1],obj.sz);
                    else
                        % A(r,c,f) with each of r,c,f either ':' or a contiguous range.
                        [st,ct] = StackH5.parse_idx(obj.sz,idx);
                        a = obj.read(st,ct);
                    end
                    % Support chained indexing like A(:,:,k)(i,j)
                    if numel(S)>1
                        a = builtin('subsref',a,S(2:end));
                    end
                otherwise
                    error('Unsupported subscript type.');
            end
        end

        function a = read(obj,st,ct)
            % Read a contiguous block from HDF5, apply flips, return double.
            % st : [startH startW startF]
            % ct : [countH countW countF]
            a = h5read(obj.f,obj.d,st,ct,[1 1 1]); % stride 1 in all dims
            if obj.flipx, a = fliplr(a); end
            if obj.flipy, a = flipud(a); end
            a = double(a);
        end

        function c = class(~)
            % Pretend to be 'double' so many MATLAB ops treat it like a normal array.
            c = 'double';
        end

        function a = double(obj)
            % Full materialization helper: read entire dataset (may be large!)
            % Prefer streaming functions below for big data.
            a = obj.read([1 1 1],obj.sz);
        end

        function m = mean(obj,dim,flag,omit)
            % Streaming mean, specialized for:
            % - mean(A,'all')
            % - mean(A,3)
            % Falls back to full materialization otherwise.
            if nargin<2, dim = 'all'; end

            if isequal(dim,'all')
                % Mean over all elements without loading entire stack at once.
                % Process in frame chunks of ~1000 by default.
                bs = [obj.sz(1) obj.sz(2) 1e3];
                % Accumulate sum and count
                [sumAll,n] = StackH5.stream_all(obj,@(x)[sum(x,'all'),numel(x)],bs);
                m = sumAll/n;
                % 'omitmissing' is accepted but data are numeric, no NaNs introduced here.
                if nargin>=4 && isequal(omit,'omitmissing'); end
                return
            end

            if isnumeric(dim) && dim==3
                % Mean across frames, returns HxW
                bs = [obj.sz(1) obj.sz(2) 1e3];
                [sumHW,c] = StackH5.stream_3(obj,@(x)sum(x,3),bs);
                m = sumHW./c;
                return
            end

            % Other cases (e.g., mean(A,1) or mean(A,2)) -> materialize.
            m = mean(double(obj),dim);
        end

        function v = std(obj,~,~,scope)
            % Streaming std for 'all' elements: std(A,[],'all')
            % Other scopes/materializations fall back to double(obj).
            if nargin<4, scope = 'all'; end
            if ~isequal(scope,'all')
                v = std(double(obj),[],scope);
                return
            end

            % One-pass for mean, second pass for sum of squares
            bs = [obj.sz(1) obj.sz(2) 1e3];
            [sumAll,n] = StackH5.stream_all(obj,@(x)[sum(x,'all'),numel(x)],bs);
            mu = sumAll/n;

            fSS = @(x) sum((x-mu).^2,'all');
            ss = StackH5.stream_all(obj,fSS,bs);
            v = sqrt(ss/(n-1));
        end

        function obj2 = fliplr(obj)
            % Return a view with horizontal flip flag toggled.
            obj2 = StackH5(obj.f,obj.d,obj.sz,~obj.flipx,obj.flipy);
        end

        function obj2 = flipud(obj)
            % Return a view with vertical flip flag toggled.
            obj2 = StackH5(obj.f,obj.d,obj.sz,obj.flipx,~obj.flipy);
        end
    end

    methods (Static, Access=private)
        function [st,ct] = parse_idx(sz,idx)
            % Convert MATLAB indices to HDF5 start/count for contiguous ranges.
            a = StackH5.span_to_st_ct(sz(1),idx,1);
            b = StackH5.span_to_st_ct(sz(2),idx,2);
            c = StackH5.span_to_st_ct(sz(3),idx,3);
            st = [a(1) b(1) c(1)];
            ct = [a(2) b(2) c(2)];
        end

        function s = norm_idx(n,ix)
            % Normalize ':' to full 1:n, otherwise pass through.
            if isequal(ix,':'), s = 1:n; return, end
            s = ix;
        end

        function out = span_to_st_ct(n,idx,k)
            % Map a 1-D index (possibly ':') to [start,count] for HDF5.
            v = StackH5.norm_idx(n,idx{k});
            if islogical(v), v = find(v); end
            if isempty(v), error('Empty indexing not supported.'); end
            if numel(v)==1
                out = [v 1]; return
            end
            if all(diff(v)==1) % contiguous
                out = [v(1) v(end)-v(1)+1];
                return
            end
            error('Non-contiguous indexing not supported; use ranges or loop.');
        end

        function varargout = stream_3(obj,f,bs)
            % Stream along the 3rd dim (frames) in chunks of size bs(3).
            % f(x) should reduce along 3rd dim (e.g., sum(x,3)).
            H = obj.sz(1); W = obj.sz(2); F = obj.sz(3);
            K = bs(3);
            s3 = 1:K:F; s3(end+1) = F+1;

            acc = []; cnt = [];
            for i=1:numel(s3)-1
                st = [1 1 s3(i)];
                ct = [H W s3(i+1)-s3(i)];
                x = obj.read(st,ct);
                y = f(x);
                if i==1
                    acc = y;
                    if nargout>1, cnt = ct(3); end
                else
                    acc = acc + y;
                    if nargout>1, cnt = cnt + ct(3); end
                end
            end
            varargout{1} = acc;
            if nargout>1, varargout{2} = cnt; end
        end

        function varargout = stream_all(obj,f,bs)
            % Stream over the whole volume in frame chunks, calling f(x)
            % and accumulating the result (scalar or [sum count] pair).
            H = obj.sz(1); W = obj.sz(2); F = obj.sz(3);
            K = bs(3);
            s3 = 1:K:F; s3(end+1) = F+1;

            acc = 0; cnt = 0;
            for i=1:numel(s3)-1
                st = [1 1 s3(i)];
                ct = [H W s3(i+1)-s3(i)];
                x = obj.read(st,ct);
                y = f(x);
                if numel(y)==2
                    % expecting [partialSum, partialCount]
                    acc = acc + y(1);
                    cnt = cnt + y(2);
                else
                    acc = acc + y; % scalar
                end
            end
            varargout{1} = acc;
            if nargout>1, varargout{2} = cnt; end
        end
    end
end
