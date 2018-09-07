%%
clear
clc
close all;
%% get path to SRCalibration

path2SRCal = 'E:\Data\Leuven Data\2018\03-Mar\22\PlaneCaibration';
path2Cal  ='E:\Data\Leuven Data\2018\03-Mar\22\PlaneCaibration\BeadsCalibrationZStack_1';

%% Initialize a zCalibration Object

calib = Core.Calib2D(path2Cal);

testSRCal = Core.SRCalibration(path2SRCal,calib.getCal);

%% get zCalibrationMovie

testSRCal.retrieveSRCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 80;

trackParam.euDistPx = 1; 
trackParam.ellip = 5;

testSRCal.retrieveSRCalData(detectParam,trackParam);

% calibrate
%% calc translation
refPlane = 3;
testSRCal.corrTranslation(refPlane);

testSRCal.checkAccuracy(refPlane);

%% calc rotation
testSRCal.corrRotation(refPlane);

testSRCal.checkAccuracy(refPlane);