% Add the WaveLab and WavDen functions to the MATLAB path.
addpath('../WavDen');
addpath(genpath('../Wavelab850'));

% Load the data.
datadir = '../../results/temp';
load(strcat(datadir,'/ml_in.txt'));

% Run BlockJS ("blockwise James-Stein") smoothing method from WavLab.
est = recblockJS('Augment',ml_in);

% Write the estimate to file, and quit.
csvwrite(strcat(datadir,'/ml_out.csv'),est);
exit
