function traces = checkTracesMan(traces)
% prompts to check traces manually
% myActivityTraces = checkTracesMan(myActivityTraces)
arguments
    traces ActivityTraces
end

% initialize
thistraces = traces.dFoverF;
idx_good = traces.goodNeuron_IDs;

to_del = [];
prompt = 'Next bad unit: [''k'' to stop]';

for i_trial = 1:traces.ntrials
    figure; imagesc(thistraces(:,idx_good,i_trial)); title(num2str(i_trial)); xlabel('units')
    
    
    while true
        x = input(prompt,'s');
        y = str2num(x);
        if isempty(y)
            if strcmp(x,'k'); break
            else disp('Input invalid.'); continue;
            end
        end
        
        % eliminating bad roi
        idx_good(y) = [];
    end

    
    close(gcf)
end

% store to traces object
traces.goodNeuron_IDs = idx_good;
traces.N = numel(idx_good);

end