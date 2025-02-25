function batchQuantifySNR(data, names)
if isstruct(data); datain = data.traces; fs = data.meta.framerate; period = data.T; names = data.trial_num;
else; datain = data; fs = 7.8125; period = 1440-4; end

for i_trial = 1:size(datain,2)
    thistrialtraces = reshape(datain(:,i_trial),period, size(datain,1)/period);
    
    RF_NT_NoiseTest1(thistrialtraces');
    
    FileOut = ['trial', num2str(names(i_trial))];
    savefig(FileOut)
end



end