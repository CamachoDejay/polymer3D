%%
clear
clc
close all;
%% get path to zCalibration


file.path = 'D:\Documents\Unif\PhD\2019-Data\04 - Apr\Interleaved\ZCal - CS';
file.ext  = '.ome.tif';
path2Cal  = 'D:\Documents\Unif\PhD\2019-Data\04 - Apr\Interleaved\2DCal';


detectParam.delta = 8;
detectParam.chi2 = 60;
fitZParam.deg = 6;
fitZParam.ellipRangeCal = [0.5 2]; %for calibration
fitZParam.ellipRange = [0.65 1.54];%To be used for data (we do not want to use too large values==> edge planes)
trackParam.euDistPx = 3; 
trackParam.commonPlanes = 2; %1 for extended, 2 for interleaved

%% Initialize a zCalibration Object
info.type = 'normal';
info.runMethod = 'run';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
info.zMethod    = 'Intensity';%or 'PSFE'

zCal = Core.ZCalibration(file,path2Cal,info);

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


