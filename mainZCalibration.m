%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'E:\Data\Leuven Data\2019\01 - Jan\8-01-zcalibration\zCal40nm';
path2Cal  = 'E:\Data\Leuven Data\2019\01 - Jan\8-01-zcalibration\planeCal40nm\FluoBeads40nm_488Exc_planecal_2';

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
calib = Core.MPCalibration(path2Cal,info);
calib.calc(4);

testZCal = Core.ZCalibration(path2zCal,calib.getCal,info);

%% get zCalibrationMovie

testZCal.retrieveZCalMov;

%% extract zData

detectParam.delta = 6;
detectParam.chi2 = 80;
fitZParam.deg = 6;
fitZParam.ellipRange = [0.5 2];

trackParam.euDistPx = 4; 
trackParam.commonPlanes = 1;

testZCal.retrieveZCalData(detectParam, fitZParam,trackParam);

%% ZCalibration

testZCal.zCalibrate;


%% show calib
method = 'spline';
testZCal.showZCalibration(method);

%% test Calibration
%fittingType = 'poly';
fittingType = 'spline';
testZCal.evalAccuracy(fittingType);
%% Save cal

testZCal.save;


