function [y, new_length] = downsampletraces(data, x)

if isstruct(data); traces = traceFormat(data.tracesdn, data.Ldn);
else traces = data;
end

downsample()
new_length = size(y,1);


end