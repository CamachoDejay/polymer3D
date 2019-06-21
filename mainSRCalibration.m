%%
clear
clc
close all;
%% get path to SRCalibration

file.path = 'E:\leuven data\05 May\13\2D Cal';
file.ext  = '.ome.tif';
path2Cal  = 'E:\leuven data\05 May\13\2D Cal';


%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'run';
info.frame2Load = 'all';
info.fitMethod = 'Phasor';
info.zMethod   = 'PSFE';

testSRCal = Core.SRCalibration(file,path2Cal,info);

%% get zCalibrationMovie

testSRCal.retrieveSRCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 40;
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