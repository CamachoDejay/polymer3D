%%
clear
clc
close all;
%% get path to zCalibration

path2zCal = 'E:\Data\Leuven Data\2019\04 - April\3\ZCal - CS';
path2Cal  = 'E:\Data\Leuven Data\2019\04 - April\3\2DCal';

detectParam.delta = 6;
detectParam.chi2 = 60;
fitZParam.deg = 6;
fitZParam.ellipRangeCal = [0.5 2]; %for calibration
fitZParam.ellipRange = [0.65 1.54];%To be used for data (we do not want to use too large values==> edge planes)
trackParam.euDistPx = 3; 
trackParam.commonPlanes = 2; %1 for extended, 2 for interleaved

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);
calib.save;
zCal = Core.ZCalibration(path2zCal,path2Cal,info);

%% get zCalibrationMovie

zCal.retrieveMovies;

%% extract zData
zCal.retrieveZCalData(detectParam, fitZParam,trackParam);

%% ZCalibration

zCal.zCalibrate;


%% show calib
method = 'spline';

zCal.showZCalibration(method);

%% test Calibration
%fittingType = 'poly';
fittingType = 'spline';
zCal.evalAccuracy(fittingType);
%% Save cal

zCal.save;


