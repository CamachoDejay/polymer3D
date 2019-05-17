clc 
clear 
close all;


path2ZCal = 'F:\Data\Leuven Data\2019\05 - May\Viscosity\Calibration';
path2SRCal = 'F:\Data\Leuven Data\2019\05 - May\Viscosity\Calibration';

path2File = 'F:\Data\Leuven Data\2019\05 - May\Viscosity\test';
path2Cal  =  'F:\Data\Leuven Data\2019\05 - May\Viscosity\Calibration';


detectParam.delta = 6;
detectParam.chi2 = 60;

%% MP Cal
info.type = 'normal';
info.runMethod = 'run';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
info.zMethod = 'PSFE';
% calib = Core.MPPlaneCalibration(path2Cal,info);
% calib.retrieveMovies;
% calib.calcIndivCal;
% calib.calcCombinedCal;
% calib.showCal(1);
% calib.save;

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


