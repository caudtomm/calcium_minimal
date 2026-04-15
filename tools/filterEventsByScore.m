
function events = filterEventsByScore(events,v,quant,method)
arguments
    events cell
    v ExperimentViewer
    quant double
    method = 'linear'
end

bda    = BaselineDriftAnalysis(v);
scores = bda.getContributionScores(method);
nSubjects = numel(scores);

idx = cell(nSubjects,1);
for sj = 1:nSubjects
    dt = scores{sj};
    
    min_val = quantile(dt, quant(1));
    max_val = quantile(dt, quant(2));
    
    idx{sj} = dt<=max_val & dt>=min_val;
end