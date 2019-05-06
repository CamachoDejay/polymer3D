clc 
clear 
close all;


path2ZCal = 'E:\Data\Leuven Data\2019\04 - April\3\ZCal - 50';
path2SRCal = 'E:\Data\Leuven Data\2019\04 - April\3\2DCal';

path2File = 'E:\Data\Leuven Data\2019\04 - April\3\Spirals';
path2Cal = 'E:\Data\Leuven Data\2019\04 - April\3\2DCal';

detectParam.delta = 6;
detectParam.chi2 = 60;

%% MP Cal
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);

%% create experiments

trackingExp = Core.TrackingExperiment(path2File,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% get TrackingData
trackParam.euDistXY = 2000;
trackParam.euDistZ  = 2000;
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam,val2Use);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% save Data

trackingExp.saveData;


