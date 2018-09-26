%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\ZCal - NormObjCorr';
path2Cal  = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\2DCal - normObjCorrPSFE\zStackFluoBeads200_S3_270618__1';

%% Initialize a zCalibration Object

calib = Core.MPCalibration(path2Cal);

testZCal = Core.ZCalibration(path2zCal,calib.getCal);

%% get zCalibrationMovie

testZCal.retrieveZCalMov;

%% extract zData

detectParam.delta = 6;
detectParam.chi2 = 80;
fitZParam.deg = 6;
fitZParam.ellipRange = [0.5 2];

trackParam.euDistPx = 6; 
trackParam.commonPlanes = 2;

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


