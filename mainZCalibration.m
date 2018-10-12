%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'E:\Data\Leuven Data\2018\ZHao\181010 - calibration\ZCalibration\PSFE';
path2Cal  = 'E:\Data\Leuven Data\2018\ZHao\181010 - calibration\2DCal\200nmFluoBeadsCal_1';

%% Initialize a zCalibration Object
info.type = 'transmission';
calib = Core.MPCalibration(path2Cal);

testZCal = Core.ZCalibration(path2zCal,calib.getCal,info);

%% get zCalibrationMovie

testZCal.retrieveZCalMov;

%% extract zData

detectParam.delta = 6;
detectParam.chi2 = 40;
fitZParam.deg = 3;
fitZParam.ellipRange = [0.7 1.42];

trackParam.euDistPx = 6; 
trackParam.commonPlanes = 2;

testZCal.retrieveZCalData(detectParam, fitZParam,trackParam);

%% ZCalibration

testZCal.zCalibrate;


%% show calib
method = 'poly';
testZCal.showZCalibration(method);

%% test Calibration
fittingType = 'poly';
%fittingType = 'spline';
testZCal.evalAccuracy(fittingType);
%% Save cal

testZCal.save;


