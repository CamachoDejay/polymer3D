clc
clear
close all
%% FilePath
%folder where folder of file are contained
file.path = 'C:\Users\MIPP\Desktop\2DCal\28-01-2019';
file.ext ='.ome.tif';

%% Initialize a zCalibration Object
info.type = 'normal';
info.nChan = 4;
info.runMethod = 'load';
info.frame2Load = 'all';

MPCal = Core.MPPlaneCalibration(file,info);

%% retrieve Cal Movie

MPCal.retrieveMovies;

%% calculate Cal
MPCal.calcIndivCal;

%% Combine Calibrations into 1
 
MPCal.calcCombinedCal;

%% test

cal = MPCal.getCal;

%% show cal
MPCal.showCal(1);

%%
camConfig = MPCal.determineCAMConfig;


%% 

MPCal.offTarget;