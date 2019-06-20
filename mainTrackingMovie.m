%%
clc;
close all;
clear;

path2ZCal = [];
path2SRCal = [];

file.path = 'F:\Data\Leuven Data\2019\Johannes\TrackingHis\Sample A_50nm_1000000_5K_1mg\1';
file.ext  = '.his';
path2Cal  =[];

detectParam.delta = 6;
detectParam.chi2 = 20;
roiRad = 7;
info.runMethod = 'load';
info.type = 'normal';
info.zMethod = 'PSFE';
info.fitMethod = 'Phasor';
%%  
MPTrackMov = Core.MPTrackingMovie(file,path2Cal,info,path2SRCal,path2ZCal);
MPTrackMov.calibrate;
MPTrackMov.giveInfo

%% test detection

 MPTrackMov.findCandidatePos(detectParam,1);
 MPTrackMov.showCandidate(1);

%% Detection
%find candidate
MPTrackMov.findCandidatePos(detectParam);

%fit position
MPTrackMov.SRLocalizeCandidate(roiRad);
%% Data correction
rot = true;
refPlane = 4;
MPTrackMov.applySRCal(rot,refPlane);

%% e-Z transformation
MPTrackMov.applyZCal;

%% Plane consolidation

MPTrackMov.consolidatePlanes


%% Super resolve
val2Use = 'bestFocus';
MPTrackMov.superResolve();
%% plot
frames = 1:10;

MPTrackMov.showCorrLoc();

%% showFrame

%MPTrackMov.showFrame(80,5);
%MPTrackMov.showParticle;

%% tracking
trackParam.radius = 1000;
trackParam.memory = 3;
MPTrackMov.trackParticle(trackParam);
traces = MPTrackMov.getTraces;
%% plot
MPTrackMov.showTraces;

%% eval accuracy
MPTrackMov.evalAccuracy('x');

%% get RMSD

[RMSD,D] = MPTrackMov.getRMSD('3D');

%% Intensity and SNR
[int,SNR] = MPTrackMov.getAvgIntensity;

%% getTraces 3D
frame = 1:87;
idx2Trace = 2;
ROIradius = 12;
frameRate = 10;
scaleBar  = 1; %in nm 
MPTrackMov.getTracesMovie(frame,idx2Trace,ROIradius,frameRate,scaleBar);
%MPTrackMov.getTraces3DMovie(frame,idx2Trace,ROIradius,frameRate);
%MPTrackMov.getPartMovie(frame,idx2Trace,ROIradius,frameRate);
