%%
clear
clc
close all;
%% get path to SRCalibration

path2SRCal = 'D:\Documents\Unif\PhD\2019-Data\04 - Apr\extended\test\2DCal';
path2Cal  = 'D:\Documents\Unif\PhD\2019-Data\04 - Apr\extended\test\2DCal';


%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod = 'Phasor';
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);

testSRCal = Core.SRCalibration(path2SRCal,calib.getCal,info);

%% get zCalibrationMovie

testSRCal.retrieveSRCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 60;
detectParam.consThresh = 3;

trackParam.commonPlanes = 1; 
trackParam.euDistPx = 3;

testSRCal.retrieveSRCalData(detectParam,trackParam);

% calibrate
%% calc translation
refPlane = 4;
testSRCal.corrTranslation(refPlane);

testSRCal.checkAccuracy(refPlane);

%% calc rotation
testSRCal.corrRotation(refPlane);

testSRCal.checkAccuracy(refPlane);