clc 
clear 
close all;


path2ZCal  = 'P:\Belgium\20190514\3D';
path2SRCal = 'P:\Belgium\20190514\3D';

path2File  = 'P:\Belgium\20190514\0%Glycerol\10 ms';
path2Cal   =  'P:\Belgium\20190514\2D';


detectParam.delta = 6;
detectParam.chi2 = 60;

%% MP Cal
info.type = 'normal';
info.runMethod = 'load';
info.frame2Load = 'all';
info.fitMethod  = 'Phasor';
info.zMethod = 'PSFE';
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);
calib.save;

%% create experiments

trackingExp = Core.TrackingExperiment(path2File,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% get TrackingData
trackParam.radius  = 1000;
trackParam.memory  = 3;
val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam,val2Use);

traces = trackingExp.getTraces3D;


%% Get Intensity

%[int,SNR] = trackingExp.getAvgIntensity;

%% save Data

trackingExp.saveData;

fileName = [path2File filesep 'fullExperiment.mat'];
save(fileName,'trackingExp','-v7.3');


