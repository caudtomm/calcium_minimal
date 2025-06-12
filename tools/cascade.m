function [spike_probs, noise_lvls, spikes] = cascade(data,envpath,framerate,pl)
    arguments
        data % numerical [t, cells, trials]
        envpath = 'C:\Users\caudtomm\AppData\Local\miniconda3\envs\Cascade\'
        framerate = 7.666
        pl = false
    end
    % prepare python environment
    terminate(pyenv)
    pyenv('Version', [envpath 'python'], 'ExecutionMode','OutOfProcess')
    setenv('PATH', [getenv('PATH') ...
    ';' envpath 'Library\bin' ...
    ';' envpath 'Library\lib']);
    pyrun('import os')
    pyrun('import numpy as np')
    pyrun('import matplotlib.pyplot as plt')
    pyrun('import glob')
    pyrun('import scipy.io as sio')
    pyrun('import glob')
    pyrun('import ruamel.yaml as yaml')
    yaml = pyrun('yaml.YAML(typ=''rt'')');
    pyrun('from cascade2p import cascade')
    pyrun('from cascade2p.utils import plot_dFF_traces, plot_noise_level_distribution, plot_noise_matched_ground_truth')
    pyrun('from cascade2p import checks')
    pyrun('checks.check_packages()')

    % convert data to numpy array
    dataNdArray = py.numpy.array(data);

    noise_levels = pyrun('plot_noise_level_distribution(traces,frame_rate)');
end