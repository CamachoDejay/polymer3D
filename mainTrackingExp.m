clc 
clear 
close all;


path2ZCal  = [];
path2SRCal = [];

path2File  = 'M:\leuven data\05 May\17 camera at different plane fluo\the 2 cameras with different plane\3W';
path2Cal   =  'M:\leuven data\05 May\17 camera at different plane fluo\the 2 cameras with different plane\2D';


detectParam.delta = 6;
detectParam.chi2 = 120;

%% MP Cal
info.type = 'normal';
info.runMethod = 'run';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
info.zMethod = 'Intensity';
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);
calib.save;

%% create experiments

trackingExp = Core.TrackingExperiment(path2File,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% get TrackingData
trackParam.radius  = 1000;
trackParam.memory  = 3;
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam,val2Use);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% save Data

trackingExp.saveData;


