%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'D:\Documents\Unif\PhD\2019-Data\01-Jan\29\ZCal';
path2Cal  = 'D:\Documents\Unif\PhD\2019-Data\01-Jan\29\PlaneCal';

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);
zCal = Core.ZCalibration(path2zCal,calib.getCal,info);

%% get zCalibrationMovie

zCal.retrieveMovies;

%% extract zData

detectParam.delta = 6;
detectParam.chi2 = 60;
fitZParam.deg = 6;
fitZParam.ellipRangeCal = [0.5 2]; %for calibration
fitZParam.ellipRange = [0.625 1.6];%To be used for data (we do not want to use too large values==> edge planes)
trackParam.euDistPx = 4; 
trackParam.commonPlanes = 1;

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


