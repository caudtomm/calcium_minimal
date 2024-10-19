function traceout = traceFormat(tracein,L)
% transforms input traces from 2D to 3D format or vice versa.
% traceout = traceFormat(tracein,L)
% 
% tracein : 1-3 D input trace matrix
% (L) : length of one trial in frames. size(tracein,1) must be a multiple
% of L




switch sum(size(tracein)>1)
    case 1
        tracein = tracein(:);
        mode = 'to3d';
    case 2
        mode = 'to3d';
    case 3
        mode = 'to2d';
    otherwise
        error('Too many dimensions!')
end

traceout = [];
switch mode
    case 'to2d'
        for i = 1:size(tracein,3)
            tmp = tracein(:,:,i);
            traceout(:,i) = tmp(:);
        end
        
    case 'to3d'
        if ~exist('L','var') || isempty(L); L = size(tracein,1); end
        n = size(tracein,1)/L;
        
        for i = 1:size(tracein,2)
            tmp = tracein(:,i);
            traceout(:,:,i) = reshape(tmp, L, n);
        end
end
        




end