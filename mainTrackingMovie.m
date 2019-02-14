%%
clc;
close all;
clear;

path2ZCal = 'D:\Documents\Unif\PhD\2019-Data\01-Jan\29\ZCal';
path2SRCal = 'D:\Documents\Unif\PhD\2019-Data\01-Jan\29\ZCal';
path2Cal = 'D:\Documents\Unif\PhD\2019-Data\01-Jan\29\PlaneCal';

path2File = 'D:\Documents\Unif\PhD\2019-Data\01-Jan\29\XYZCAL\FluoBeads200nm_488Exc_XY_3';

detectParam.delta = 6;
detectParam.chi2 = 80;
info.runMethod = 'load';
info.type = 'normal';
%%
calib = Core.MPPlaneCalibration(path2Cal,info);
calib.retrieveMovies;
calib.calcIndivCal;
calib.calcCombinedCal;
calib.showCal(1);
%%

MPTrackMov = Core.MPTrackingMovie(path2File,calib.getCal,info,path2SRCal,path2ZCal);
MPTrackMov.calibrate;
%% Detection

MPTrackMov.giveInfo
%find candidate
MPTrackMov.findCandidatePos(detectParam);

%fit position
MPTrackMov.SRLocalizeCandidate;

%% Data correction
rot = true;
refPlane = 3;
MPTrackMov.applySRCal(rot,refPlane);
%% e-Z transformation
MPTrackMov.applyZCal;

%% Plane consolidation

MPTrackMov.consolidatePlanes
%% Super resolve
val2Use = 'bestFocus';
MPTrackMov.superResolve(val2Use);
%% plot
frames = 1:500;

MPTrackMov.showCorrLoc();

%% showFrame

%MPTrackMov.showFrame(80,5);
%MPTrackMov.showParticle;

%% tracking
trackParam.euDistXY = 800;
trackParam.euDistZ  = 400;
MPTrackMov.trackParticle(trackParam);
traces = MPTrackMov.getTraces;
%% plot
MPTrackMov.showTraces;

%% eval accuracy
MPTrackMov.evalAccuracy

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
