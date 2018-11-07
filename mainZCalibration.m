%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'E:\Data\Leuven Data\2018\10-Oct\23\ZCal\200nmBeads';
path2Cal  = 'E:\Data\Leuven Data\2018\10-Oct\23\2DCal\FluoBeads200nm_1';

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
calib = Core.MPCalibration(path2Cal);
calib.calc(1);

testZCal = Core.ZCalibration(path2zCal,calib.getCal,info);

%% get zCalibrationMovie

testZCal.retrieveZCalMov;

%% extract zData

detectParam.delta = 6;
detectParam.chi2 = 60;
fitZParam.deg = 4;
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


