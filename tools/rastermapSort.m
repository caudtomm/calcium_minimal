 function [data_sorted, sortIdx, hf] = rastermapSort(data,envpath, locality, time_lag_w, pl)
    arguments
        data % numerical [t, cells, trials]
        envpath = 'C:\Users\caudtomm\AppData\Local\miniconda3\envs\rastermap\python'
        locality = .1
        time_lag_w = 50
        pl = false
    end
    % prepare python environment
    terminate(pyenv)
    pyenv('Version', envpath, 'ExecutionMode','OutOfProcess')
    pyrun("from rastermap import Rastermap")

    % prepare data
    data_input = [];
    for i = 1:size(data,3); data_input = [data_input; data(:,:,i)]; end
    data_input = nanzscore(data_input)';
    idx = isnan(data_input(1,:));
    data_input(:,idx) = [];

    % convert data to numpy array
    dataNdArray = py.numpy.array(data_input);

    % create Rastermap model
    rmModel = pyrun(["model = Rastermap(locality=",num2str(locality),...
        ", time_lag_window=",num2str(time_lag_w),").fit(spks)"], ...
        "model", spks=dataNdArray);

    % convert indices back to MATLAB array
    sortIdx = int16(py.memoryview(rmModel.isort.data)) + 1;

    % sort original data
    data_sorted = data(:,sortIdx,:);

    % optional plotting
    hf = [];
    if ~pl; return; end
    
    data_input = [];
    for i = 1:size(data,3); data_input = [data_input; data(:,:,i)]; end
    data_input = nanzscore(data_input)';

    data_output = data_input(sortIdx,:);

    hf = figure;
    subplot(211); imagesc(data_input); colormap('gray')
    title('input data'); ylabel('neuron #')
    subplot(212); imagesc(data_output); colormap('gray')
    title('sorted data'); ylabel('neuron #'); xlabel('frames')
end