clc
clear
close all
%% FilePath
%folder where folder of file are contained
path2MPCal = 'E:\Data\Leuven Data\2019\01 - Jan\8-01-zcalibration\planeCal200nm';

%% Initialize a zCalibration Object
info.type = 'normal';
info.nChan = 4;
info.runMethod = 'load';
info.frame2Load = 'all';

MPCal = Core.PlaneCalibration(path2MPCal,info);

%% retrieve Cal Movie

MPCal.retrieveMPCalibration;

%% calculate Cal
MPCal.calcIndivCal;

%% Combine Calibrations into 1
 
MPCal.calcCombinedCal;