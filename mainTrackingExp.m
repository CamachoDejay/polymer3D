clc 
clear 
close all;


path2ZCal = 'F:\Data\Leuven Data\2019\04 - April\4\ZCal CS';
path2SRCal = 'F:\Data\Leuven Data\2019\04 - April\4\2DCal';

path2File = 'F:\Data\Leuven Data\2019\04 - April\4\SPiral';
path2Cal = 'F:\Data\Leuven Data\2019\04 - April\4\2DCal';

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
trackParam.euDistXY = 1500;
trackParam.euDistZ  = 300;
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam,val2Use);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% save Data

trackingExp.saveData;


