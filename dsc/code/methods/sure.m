% Add the WaveLab and WavDen functions to the MATLAB path.
addpath('../WavDen');
addpath(genpath('../Wavelab850'));

% Load the data.
datadir = '../../results/temp';
load(strcat(datadir,'/ml_in.txt'));

% Run "sureshrink" (Stein's unbiased risk estimation) smoothing method
% from WavLab.
est = recsure(ml_in);

% Write the estimate to file, and quit.
csvwrite(strcat(datadir,'/ml_out.csv'),est);
exit
