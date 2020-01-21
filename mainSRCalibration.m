%%
clear
clc
close all;
%% get path to SRCalibration

file.path = 'D:\Dropbox\4Apr Data\2DCal';
file.ext  = '.ome.tif';
path2Cal  = 'D:\Dropbox\4Apr Data\2DCal';


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
detectParam.chi2 = 60;
detectParam.consThresh = 4;

trackParam.commonPlanes = 1; 
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