clc 
clear 
close all;


path2ZCal = 'E:\Data\Leuven Data\2019\04 - April\3\ZCal - CS';
path2SRCal = 'E:\Data\Leuven Data\2019\04 - April\3\2DCal';

path2File = 'E:\Data\Leuven Data\2018\06-June\29\5K - 0_25mgmL\';
path2Cal = 'E:\Data\Leuven Data\2019\04 - April\3\2DCal';

detectParam.delta = 6;
detectParam.chi2 = 60;

%% MP Cal

calib = Core.MPCalibration(path2Cal);

%% create experiments

trackingExp = Core.TrackingExperiment(path2File,calib.getCal,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% get TrackingData

trackingExp.retriveTrackData;