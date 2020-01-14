clc 
clear 
close all;

path2ZCal = [];
path2SRCal = 'F:\Boris\2019\ViscosityData\20190627\2D Cal\200 nm';

file.path  = 'F:\Boris\2019\ViscosityData\20190627\500nm-Fluo';
file.ext   = '.ome.tif';
path2Cal = 'F:\Boris\2019\ViscosityData\20190627\2D Cal\200 nm';
dimension = '3D';
detectParam.delta = 6;
detectParam.chi2 = 80;

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
frame = 20;
testMov = trackingExp.trackMovies.mov1;
testMov.findCandidatePos(detectParam,frame);
testMov.showCandidate(frame);

%% get TrackingData
trackParam.radius  = 1000;
trackParam.memory  = 3;

val2Use = 'bestFocus';
trackingExp.retrieveTrackData(detectParam,trackParam);

traces = trackingExp.getTraces3D;


%% Get Intensity

[int,SNR] = trackingExp.getAvgIntensity;

%% Get MSD

[MSD,~] = trackingExp.getMSD(dimension);
%% show traces

trackingExp.showTraces(1);
%% save Data

trackingExp.saveData;


