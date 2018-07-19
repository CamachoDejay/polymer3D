%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'E:\Data\Leuven Data\2018\03-Mar\22\ZCalibration';
path2Cal  ='E:\Data\Leuven Data\2018\03-Mar\22\PlaneCaibration\BeadsCalibrationZStack_1';

%% Initialize a zCalibration Object

calib = Core.calib2D(path2Cal);

testZCal = Core.zCalibration(path2zCal,calib.getCal);

%% get zCalibrationMovie

testZCal.retrieveZCalMov;

%% extract zData
detectParam.delta = 6;
detectParam.chi2 = 80;
testZCal.retrieveZCalData(detectParam)

%% ZCalibration

testZCal.zCalibrate;


%% show calib

testZCal.showZCalibration