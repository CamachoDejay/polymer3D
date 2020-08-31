clc 
clear 
close all;
%calibration info
path2ZCal = [];
path2SRCal = 'D:\Documents\Unif\PhD\2020-Data\08 - August\Johannes\Elastic polymer\2DCal';

%file info
file.path  = 'D:\Documents\Unif\PhD\2020-Data\08 - August\Johannes\Elastic polymer\Mov3 - Ok';
file.ext   = '.ome.tif';
path2Cal = 'D:\Documents\Unif\PhD\2020-Data\08 - August\Johannes\Elastic polymer\2DCal';
dimension = '3D';

%detection parameter
detectParam.delta = 6;
detectParam.chi2  = 80;
detectParam.consThresh = 4;
%tracking parameter
trackParam.radius  = 1000;%nm
trackParam.memory  = 3;

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'run'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.calibrate = false; %true to recalibrate;

%% create experiments

trackingExp = Core.TrackingExperiment(file,path2Cal,info,path2SRCal,path2ZCal);

%% get Movies

trackingExp.retrieveMovies;

%% test detection parameters
frame =99;
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam,frame);
testMov.showCandidate(frame);

%% get TrackingData


val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% Get MSD

%[MSD,~] = trackingExp.getMSD(dimension);
%% show traces

trackingExp.showTraces(1);
%% save Data

trackingExp.saveData;


