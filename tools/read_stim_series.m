function [stim_series, ch2stimtype, ch2stimtype_map] = read_stim_series(filename)

% filename = 'C:\Users\caudtomm\Desktop\2021_invivo_data\TC_210505_TC0004_01_sxaDp_odorexp001_RPB301AG__00001'

M = importdata(filename);
ntrials = floor((numel(M)-1)/4);
trial_lines = 2:4:numel(M); if length(trial_lines)>ntrials; trial_lines = trial_lines(1:ntrials); end
stim_series = [];

for i = 1:length(trial_lines)
    trial_num = i;
    a = strsplit(M{trial_lines(i)},'\t');
    channel_id = str2num(a{3});
    stim_on_fr = floor(30*7.67);
    stim_off_fr = floor(50*7.67);
    stim_series = [stim_series; trial_num, channel_id, stim_on_fr, stim_off_fr];
end
clear M

FileIn = strcat(filename, '_stimtypes.csv');
if ~exist(FileIn,'file')
    error('no stimtype map file found!')
end

N = importdata(FileIn);
ch2stimtype_map = {};
for i = 2:numel(N)
    x = strsplit(N{i},',');
    ch2stimtype_map{end+1} = x{2};
end
ch2stimtype_map = ch2stimtype_map';

ch2stimtype = {};
for i = 1:ntrials
    ch2stimtype{end+1} = ch2stimtype_map{stim_series(i,2)};
end
ch2stimtype = ch2stimtype';