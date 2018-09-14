%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'E:\Data\Leuven Data\2018\03-Mar\22\ZCalibration';
path2Cal  = 'E:\Data\Leuven Data\2018\03-Mar\22\PlaneCaibration\BeadsCalibrationZStack_1';

%% Initialize a zCalibration Object

calib = Core.MPCalibration(path2Cal);

testZCal = Core.ZCalibration(path2zCal,calib.getCal);

%% get zCalibrationMovie

testZCal.retrieveZCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 80;
fitZParam.deg = 6;
fitZParam.ellipRange = [0.6 2];

trackParam.euDistPx = 2; 
trackParam.commonPlanes = 2;

testZCal.retrieveZCalData(detectParam, fitZParam,trackParam);

%% ZCalibration

testZCal.zCalibrate;


%% show calib

testZCal.showZCalibration

%% test Calibration
fittingType = 'poly';
%fittingType = 'spline';
testZCal.evalAccuracy(fittingType);
%% Save cal

testZCal.save;


