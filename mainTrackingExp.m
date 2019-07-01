clc 
clear 
close all;

path2ZCal = 'F:\Data\Leuven Data\2019\04 - April\4\ZCal CS';
path2SRCal = 'F:\Data\Leuven Data\2019\04 - April\4\2DCal';

file.path  = 'F:\Data\Leuven Data\2019\04 - April\4\XYZ - OD50\Z';
file.ext   = '.ome.tif';
path2Cal = 'F:\Data\Leuven Data\2019\04 - April\4\2DCal';
dimension = '3D';
detectParam.delta = 6;
detectParam.chi2 = 24;

%% Storing info about the file
info.type = 'normal';
info.runMethod = 'run';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
info.zMethod = 'PSFE';

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


