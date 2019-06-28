clc 
clear 
close all;


path2ZCal  = 'F:\Data\Leuven Data\2019\05 - May\Viscosity\Calibration';
path2SRCal = 'F:\Data\Leuven Data\2019\05 - May\Viscosity\Calibration';

file.path  = 'F:\Data\Leuven Data\2019\05 - May\Viscosity\glycerol';
file.ext   = '.ome.tif';
path2Cal   =  'F:\Data\Leuven Data\2019\05 - May\Viscosity\Calibration';

dimension = '3D';
detectParam.delta = 6;
detectParam.chi2 = 60;

%% Storing info about the file
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
info.zMethod = 'Intensity';

%% create experiments

trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;


%% get TrackingData
trackParam.radius  = 1000;
trackParam.memory  = 3;

val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% Get MSD

[MSD,~] = trackingExp.getRMSD(dimension);

%% save Data

trackingExp.saveData;


