classdef StackH5Multi
    properties (SetAccess=private)
        files string        % shard file paths
        lens  double        % frames per shard
        offs  double        % starting frame index of each shard (1-based)
        sz    (1,3) double  % [H W F]
        flipx logical = false
        flipy logical = false
    end
    methods
        function obj = StackH5Multi(files,lens,sz,flipx,flipy)
            obj.files = files(:);
            obj.lens  = lens(:);
            obj.offs  = cumsum([1; lens(1:end-1)]);
            obj.sz    = sz;
            if nargin>=4, obj.flipx = logical(flipx); end
            if nargin>=5, obj.flipy = logical(flipy); end
        end
        function varargout = size(obj, dim)
            sz = obj.sz;
            if nargin>1
                s = dim<=numel(sz) * sz(min(dim,end)) + (dim>numel(sz));
                varargout = {s};
            else
                if nargout<=1, varargout = {sz};
                else, out = [sz, ones(1, nargout-numel(sz))]; varargout = num2cell(out(1:nargout));
                end
            end
        end
        function n = ndims(~), n = 3; end
        function c = class(~), c = 'double'; end

        function a = subsref(obj,S)
            switch S(1).type
                case '.'
                    a = builtin('subsref',obj,S);

                case '()'
                    idx = S(1).subs;

                    % normalize colon flags
                    iscolon = @(x) (ischar(x) && strcmp(x,':')) || (isstring(x) && x==":");
                    
                    % Fast path: scalar frame index in 3rd dim
                    if numel(idx)==3 && (iscolon(idx{1}) || isscalar(idx{1})) ...
                                    && (iscolon(idx{2}) || isscalar(idx{2})) ...
                                    && isscalar(idx{3}) && ~iscolon(idx{3})
                        % Normalize dims 1,2 to full span if ':'
                        if iscolon(idx{1}), r = [1 obj.sz(1)]; else, r = [idx{1} 1]; end
                        if iscolon(idx{2}), c = [1 obj.sz(2)]; else, c = [idx{2} 1]; end
                        f = idx{3};
                        a = obj.read([r(1) c(1) f],[r(2) c(2) 1]);  % returns 2D
                    else
                        % General path
                        [st,ct] = parse_idx(obj.sz,idx);
                        a = obj.read(st,ct);
                    end

                    % Support chaining
                    if numel(S)>1
                        a = builtin('subsref',a,S(2:end));
                    end

                otherwise
                    error('Unsupported subscript type.');
            end
        end


       function a = read(obj,st,ct)
            % st = [row_start col_start frame_start], ct = [nrows ncols nframes]
            f0 = st(3);
            k  = ct(3);
            H  = ct(1);
            W  = ct(2);

            % 1) Fast path: single frame -> return 2D slice (small, no big allocs)
            if k == 1
                f = f0;
                i = find( f >= obj.offs & f <= (obj.offs + obj.lens - 1), 1, 'first' );
                if isempty(i), error('Frame index out of range.'); end
                loc = f - obj.offs(i) + 1;                                     % 1-based within shard
                blk = h5read(obj.files(i), '/stack', [st(1) st(2) loc], [ct(1) ct(2) 1]);
                if obj.flipx, blk = fliplr(blk); end
                if obj.flipy, blk = flipud(blk); end
                a = double(blk(:,:,1));                                         % 2D double
                return
            end

            % 2) Memory guard: very large request -> stream to temp HDF5, return StackH5 view
            bytes_double = double(H) * double(W) * double(k) * 8;
            memlim = get_mem_limit_bytes();  % ~1/3 of MaxPossibleArrayBytes (see helper)
            if bytes_double > memlim
                tmpf = [tempname '.h5'];
                h5create(tmpf, '/stack', [H W k], 'Datatype','uint8', ...
                    'Chunksize', [H W min(k,128)], 'Deflate', 0);

                % fill in manageable frame chunks
                K = min(512, k);                         % frames per chunk
                f_lo = f0; f_hi = f0 + k - 1;
                dst  = 1;
                while f_lo <= f_hi
                    f_end = min(f_lo + K - 1, f_hi);
                    loc_k = f_end - f_lo + 1;

                    % assemble chunk from shards
                    blk = zeros(H, W, loc_k, 'uint8');
                    for iShard = 1:numel(obj.files)
                        s0 = obj.offs(iShard); s1 = s0 + obj.lens(iShard) - 1;
                        lo = max(f_lo, s0); hi = min(f_end, s1);
                        if lo > hi, continue; end
                        loc_start = lo - s0 + 1;
                        loc_count = hi - lo + 1;
                        dst_start = (lo - f_lo) + 1;
                        sub = h5read(obj.files(iShard), '/stack', [st(1) st(2) loc_start], [H W loc_count]);
                        if obj.flipx, sub = fliplr(sub); end
                        if obj.flipy, sub = flipud(sub); end
                        blk(:,:,dst_start:dst_start+loc_count-1) = sub;
                    end

                    h5write(tmpf, '/stack', blk, [1 1 dst], [H W loc_k]);
                    dst  = dst  + loc_k;
                    f_lo = f_end + 1;
                end

                a = StackH5(tmpf, '/stack', [H W k], obj.flipx, obj.flipy);
                return
            end

            % 3) Normal multi-frame path (fits in RAM): allocate once and fill
            a  = zeros(H, W, k, 'double');
            f1 = f0 + k - 1;
            for iShard = 1:numel(obj.files)
                s0 = obj.offs(iShard);
                s1 = s0 + obj.lens(iShard) - 1;
                lo = max(f0, s0);
                hi = min(f1, s1);
                if lo>hi, continue; end
                loc_start = lo - s0 + 1;
                loc_count = hi - lo + 1;
                dst_start = lo - f0 + 1;
                blk = h5read(obj.files(iShard), '/stack', [st(1) st(2) loc_start], [H W loc_count]);
                if obj.flipx, blk = fliplr(blk); end
                if obj.flipy, blk = flipud(blk); end
                a(:,:,dst_start:dst_start+loc_count-1) = double(blk);
            end
        end

        function m = mean(obj,dim,~,omit)
            if nargin<2, dim = 'all'; end
            if isequal(dim,'all')
                bs = [obj.sz(1) obj.sz(2) 1024];
                [s,c] = stream_all(obj,@(x)[sum(x,'all'),numel(x)],bs);
                m = s/c; if nargin>=4 && isequal(omit,'omitmissing'); end
                return
            end
            if isnumeric(dim) && dim==3
                bs = [obj.sz(1) obj.sz(2) 1024];
                [s,c] = stream_3(obj,@(x)sum(x,3),bs);
                m = s./c; return
            end
            m = mean(double(obj),dim);
        end

        function v = std(obj,~,~,scope)
            if nargin<4, scope = 'all'; end
            if ~isequal(scope,'all'), v = std(double(obj),[],scope); return, end
            bs = [obj.sz(1) obj.sz(2) 1024];
            [s,c] = stream_all(obj,@(x)[sum(x,'all'),numel(x)],bs);
            mu = s/c;
            ss = stream_all(obj,@(x)sum((x-mu).^2,'all'),bs);
            v = sqrt(ss/(c-1));
        end

        function obj2 = fliplr(obj), obj2 = StackH5Multi(obj.files,obj.lens,obj.sz,~obj.flipx,obj.flipy); end
        function obj2 = flipud(obj), obj2 = StackH5Multi(obj.files,obj.lens,obj.sz,obj.flipx,~obj.flipy); end
    end
end

% ---- helpers ----
function [st,ct] = parse_idx(sz,idx)
    a = span_to_st_ct(sz(1),idx,1);
    b = span_to_st_ct(sz(2),idx,2);
    c = span_to_st_ct(sz(3),idx,3);
    st = [a(1) b(1) c(1)];
    ct = [a(2) b(2) c(2)];
    % bounds guard (using sz)
    if st(1) < 1 || st(1)+ct(1)-1 > sz(1) || ...
       st(2) < 1 || st(2)+ct(2)-1 > sz(2) || ...
       st(3) < 1 || st(3)+ct(3)-1 > sz(3)
        error('Requested index range exceeds stack size [%d %d %d].', sz(1), sz(2), sz(3));
    end
end
function out = span_to_st_ct(n,idx,k)
    v = idx{k};
    if isequal(v,':'), v = 1:n; end
    if islogical(v), v = find(v); end
    if isempty(v), error('Empty indexing not supported.'); end
    if numel(v)==1
        out = [v 1];
    elseif all(diff(v)==1)
        out = [v(1) v(end)-v(1)+1];
    else
        error('Non-contiguous indexing not supported; use ranges or loops.');
    end
end
function varargout = stream_3(obj,f,bs)
    H = obj.sz(1); W = obj.sz(2); F = obj.sz(3);
    K = bs(3); s3 = 1:K:F; s3(end+1) = F+1;
    acc = []; cnt = 0;
    for i=1:numel(s3)-1
        st=[1 1 s3(i)]; ct=[H W s3(i+1)-s3(i)];
        x = obj.read(st,ct); y=f(x);
        if i==1, acc=y; cnt=ct(3); else, acc=acc+y; cnt=cnt+ct(3); end
    end
    varargout{1}=acc; if nargout>1, varargout{2}=cnt; end
end
function varargout = stream_all(obj,f,bs)
    H = obj.sz(1); W = obj.sz(2); F = obj.sz(3);
    K = bs(3); s3 = 1:K:F; s3(end+1) = F+1;
    acc=0; cnt=0;
    for i=1:numel(s3)-1
        st=[1 1 s3(i)]; ct=[H W s3(i+1)-s3(i)];
        x = obj.read(st,ct); y=f(x);
        if numel(y)==2, acc=acc+y(1); cnt=cnt+y(2); else, acc=acc+y; end
    end
    varargout{1}=acc; if nargout>1, varargout{2}=cnt; end
end

function lim = get_mem_limit_bytes()
try
    m = memory;                        % works on Windows
    lim = max(1e9, floor(m.MaxPossibleArrayBytes/3));  % be conservative
catch
    lim = 1.5e9;                       % ~1.5 GB fallback on non-Windows
end
end
