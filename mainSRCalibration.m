%%
clear
clc
close all;
%% get path to SRCalibration

file.path = 'F:\Boris\2019\ViscosityData\20190627\2D Cal\200 nm';
file.ext  = '.ome.tif';
path2Cal  = 'F:\Boris\2019\ViscosityData\20190627\2D Cal\200 nm';


%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod = 'Phasor';
info.zMethod   = 'Intensity';

testSRCal = Core.SRCalibration(file,path2Cal,info);

%% get zCalibrationMovie

testSRCal.retrieveSRCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 80;
detectParam.consThresh = 4;

trackParam.commonPlanes = 2; 
trackParam.euDistPx = 3;

testSRCal.retrieveSRCalData(detectParam,trackParam);

% calibratewwwnBB
%% calc translation
refPlane = 4;
testSRCal.corrTranslation(refPlane);

testSRCal.checkAccuracy(refPlane);

%% calc rotation
testSRCal.corrRotation(refPlane);

testSRCal.checkAccuracy(refPlane);