%%
clear
clc
close all;
%% get path to SRCalibration

path2SRCal = 'E:\Data\Leuven Data\2019\01 - Jan\8-01-zcalibration\zCal200nm';
path2Cal  = 'E:\Data\Leuven Data\2019\01 - Jan\8-01-zcalibration\planeCal200nm\FluoBeads200nm_488Exc_planecal_2';


%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
calib = Core.MPCalibration(path2Cal,info);
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