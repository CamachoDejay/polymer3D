clc 
clear 
close all;


path2ZCal = 'E:\Data\Leuven Data\2019\04 - April\3\ZCal - CS';
path2SRCal = 'E:\Data\Leuven Data\2019\04 - April\3\2DCal';

path2File = 'E:\Data\Leuven Data\2019\04 - April\3\test\X';
path2Cal = 'E:\Data\Leuven Data\2019\04 - April\3\2DCal';

detectParam.delta = 6;
detectParam.chi2 = 60;

%% MP Cal
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);

%% create experiments

trackingExp = Core.TrackingExperiment(path2File,calib.getCal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% get TrackingData
detectParam.delta = 6;
detectParam.chi2 = 60;
trackParam.euDistXY = 800;
trackParam.euDistZ  = 400;
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam,val2Use);