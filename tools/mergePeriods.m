function periods_out = mergePeriods(periods_in,th)
% expects periods with 3 columns [trial, start, end]
% threshold period length 'th' [frames]

% initialize output
periods_out = [];

trials = unique(periods_in(:,1));
if isempty(trials); return; end

for i = 1:numel(trials)
    % isolate this trial's periods
    thistrial = trials(i);
    thisperiods = periods_in(periods_in(:,1)==thistrial,:);

    % find any periods over the threshold duration
    dur = diff(thisperiods(:,[2,3]),[],2);
    idx_long = dur >= th;

    % store only longer periods to output
    periods_out = [periods_out; thisperiods(idx_long,:)];
end

end