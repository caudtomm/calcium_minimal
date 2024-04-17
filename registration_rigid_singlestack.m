function [movie,movie_raw,stats] = registration_rigid_singlestack(FilenameRAW,FilenameCLAHE,ref_img,badperiods,framerate)

FilenameRAW = char(FilenameRAW);
FilenameCLAHE = char(FilenameCLAHE);

fspecs = getFileNameSpecs(FilenameRAW);

disp('')
disp('Linear registration:')
disp(fspecs.fname)

L = imfinfo(FilenameRAW);
meta.height = L(1).Height;
meta.width = L(1).Width;

movie = double(loadTiffStack(FilenameCLAHE));
movie_raw = double(loadTiffStack(FilenameRAW));
meta.numberframes = size(movie_raw,3);

meta.framerate = framerate;

%% Subtract bad periods
[movie, ~] = subtract_badperiods(movie,badperiods);
[movie_raw, badperiods] = subtract_badperiods(movie_raw,badperiods);

%% Subtract baseline from movie
[movie,baseline_histeq] = subtract_movie_baseline(movie,meta.framerate);
[movie_raw,baseline_raw] = subtract_movie_baseline(movie_raw,meta.framerate);

%% register
[movie,movie_raw,shift_x,shift_y] = linearAlignMovie(movie,movie_raw,ref_img);

%% add nans in the place of bad periods
movie = add_badperiods_as_nans(movie,badperiods,meta);
movie_raw = add_badperiods_as_nans(movie_raw,badperiods,meta);

%% save linear registration result

pathout = 'reg_stacks3_clahe';
if ~exist(pathout,'dir'); mkdir(pathout); end
fileout = [pathout '\' fspecs.fname '_reg.tif'];
if exist(fileout,'file'); delete(fileout); end
saveastiff(uint16(movie), fileout)  

pathout = 'reg_stacks3_raw'
if ~exist(pathout,'dir'); mkdir(pathout); end
fileout = [pathout '\' fspecs.fname '_reg.tif'];
if exist(fileout,'file'); delete(fileout); end
saveastiff(uint16(movie_raw), fileout)  

pathout = 'reg_stacks3_raw';
if ~exist(pathout,'dir'); mkdir(pathout); end
fileout = [pathout '\' fspecs.fname '_stats.m'];
if exist(fileout,'file'); delete(fileout); end
stats.shifts = [shift_x,shift_y];
stats.baseline_raw = baseline_raw;
stats.baseline_histeq = baseline_histeq;
save(fileout,'stats','-v7.3');
