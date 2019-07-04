clc 
clear 
close all;

path2ZCal = [];
path2SRCal = [];

file.path  = 'F:\Data\Leuven Data\2019\06 - June\Roger\ZPos\test 2DTrans 1000';
file.ext   = '.ome.tif';
path2Cal = 'F:\Data\Leuven Data\2019\06 - June\Roger\2Dcal\200 nm';
dimension = '3D';
detectParam.delta = 6;
detectParam.chi2 = 60;

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

%% test detection parameters
frame = 65;
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

[MSD,~] = trackingExp.getRMSD(dimension);
%% show traces

trackingExp.showTraces(1);
%% save Data

trackingExp.saveData;


