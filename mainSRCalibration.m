%%
clear
clc
close all;
%% get path to SRCalibration

path2SRCal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE';
path2Cal  ='E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE\zStackFluoBeads200_S3_270618__1';

%% Initialize a zCalibration Object

calib = Core.MPCalibration(path2Cal);

testSRCal = Core.SRCalibration(path2SRCal,calib.getCal);

%% get zCalibrationMovie

testSRCal.retrieveSRCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 80;

trackParam.commonPlanes = 2; 
trackParam.euDistPx = 5;

testSRCal.retrieveSRCalData(detectParam,trackParam);

% calibrate
%% calc translation
refPlane = 5;
testSRCal.corrTranslation(refPlane);

testSRCal.checkAccuracy(refPlane);

%% calc rotation
testSRCal.corrRotation(refPlane);

testSRCal.checkAccuracy(refPlane);