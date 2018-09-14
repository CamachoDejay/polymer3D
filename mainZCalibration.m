%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - NormObjCorr';
path2Cal  = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE\zStackFluoBeads200_S3_270618__1';

%% Initialize a zCalibration Object

calib = Core.MPCalibration(path2Cal);

testZCal = Core.ZCalibration(path2zCal,calib.getCal);

%% get zCalibrationMovie

testZCal.retrieveZCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 80;
fitZParam.deg = 6;
fitZParam.ellipRange = [0.6 1.66];

trackParam.euDistPx = 6; 
trackParam.commonPlanes = 1;

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


