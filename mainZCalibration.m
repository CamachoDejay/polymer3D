%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'D:\Documents\Unif\PhD\2019-Data\Roger\02-Feb\26\3D Cal';
path2Cal  = 'D:\Documents\Unif\PhD\2019-Data\Roger\02-Feb\26\2DCal\Crop';

detectParam.delta = 6;
detectParam.chi2 = 60;
fitZParam.deg = 6;
fitZParam.ellipRangeCal = [0.5 2]; %for calibration
fitZParam.ellipRange = [0.7 1.42];%To be used for data (we do not want to use too large values==> edge planes)
trackParam.euDistPx = 3; 
trackParam.commonPlanes = 1; %1 for extended, 2 for interleaved

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);
zCal = Core.ZCalibration(path2zCal,calib.getCal,info);

%% get zCalibrationMovie

zCal.retrieveMovies;

%% extract zData
zCal.retrieveZCalData(detectParam, fitZParam,trackParam);

%% ZCalibration

zCal.zCalibrate;


%% show calib
method = 'spline';

zCal.showZCalibration(method);

%% test Calibration
%fittingType = 'poly';
fittingType = 'spline';
zCal.evalAccuracy(fittingType);
%% Save cal

zCal.save;


