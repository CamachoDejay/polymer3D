%%
clear
clc
close all;
%% get path to SRCalibration

path2SRCal = 'E:\Data\Leuven Data\2018\10-Oct\23\2DCal';
path2Cal  ='E:\Data\Leuven Data\2018\10-Oct\23\2DCal\FluoBeads200nm_1';

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';

calib = Core.MPCalibration(path2Cal);
calib.calc(1)

testSRCal = Core.SRCalibration(path2SRCal,calib.getCal,info);

%% get zCalibrationMovie

testSRCal.retrieveSRCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 60;
detectParam.consThresh = 5;

trackParam.commonPlanes = 1; 
trackParam.euDistPx = 5;

testSRCal.retrieveSRCalData(detectParam,trackParam);

% calibrate
%% calc translation
refPlane = 2;
testSRCal.corrTranslation(refPlane);

testSRCal.checkAccuracy(refPlane);

%% calc rotation
testSRCal.corrRotation(refPlane);

testSRCal.checkAccuracy(refPlane);