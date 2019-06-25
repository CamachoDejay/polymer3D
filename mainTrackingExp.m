clc 
clear 
close all;


path2ZCal  = [];
path2SRCal = [];

file.path  = 'F:\Data\Leuven Data\2019\Johannes\TrackingHis\Sample C_50nm_1000000_1K_1mg\fastAcq';
file.ext   = '.his';
path2Cal   =  [];

dimension = '2D';
detectParam.delta = 6;
detectParam.chi2 = 24;

%% Storing info about the file
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
info.zMethod = 'Intensity';

%% create experiments

trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;


%% get TrackingData
trackParam.radius  = 1000;
trackParam.memory  = 3;

val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam,val2Use);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% Get MSD

%[MSD,~] = trackingExp.getRMSD(dimension);

%% save Data

trackingExp.saveData;


