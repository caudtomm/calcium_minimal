function [periods, len] = convertPeriods(X,reverse, len)
    arguments
        X % logical array or period matrix
        reverse logical = false;
        len double = 0
    end
    
    if ~reverse % X is logical
        X = logical(X);
        len = length(X);
    
        pstart = find(diff(X)==1)+1;
        pend = find(diff(X)==-1);

        % make sure array length is consistent
        if numel(pstart) > numel(pend)
            pend = [pend; length(X)];
        elseif numel(pstart) < numel(pend)
            pstart = [1; pstart];
        end

        % make sure periods are not staggered
        if ~isempty(pend) && ~isempty(pstart) && pend(1)<pstart(1)
            pend = [pend; length(X)];
            pstart = [1; pstart];
        end
        
        periods = [pstart, pend];
    else % X is a period matrix [:,[pstart, pend]]
        len = max([len, max(X,[],"all")]);
        periods = false(len,1);
        nps = size(X,1);
        for i = 1:nps
            thisperiod = X(i,:);
            periods(thisperiod(1):thisperiod(2)) = true;
        end
    end
end