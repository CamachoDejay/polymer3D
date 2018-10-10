%%
clc 
close all
clear 

path2ZCal = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\ZCal - NormObjCorr';
path2SRCal = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\2DCal - normObjCorrPSFE';

path2File = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\trackingCal - normObjCorr\TrackingX\FluoBeads200_TrackingX_320ms__1';
path2Cal = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\2DCal - normObjCorrPSFE\zStackFluoBeads200_S3_270618__1';

detectParam.delta = 6;
detectParam.chi2 = 80;
%%
calib = Core.MPCalibration(path2Cal);

%%

MPTrackMov = Core.MPTrackingMovie(path2File,calib.getCal,path2SRCal,path2ZCal);

%% Detection

MPTrackMov.giveInfo
%find candidate
MPTrackMov.findCandidatePos(detectParam);

%fit position
MPTrackMov.SRLocalizeCandidate;

%% Data correction
rot = true;
refPlane = 5;
MPTrackMov.applySRCal(rot,refPlane);
%% e-Z transformation
MPTrackMov.applyZCal;

%% Plane consolidation

MPTrackMov.consolidatePlanes
%% Super resolve
val2Use = 'bestFocus';
MPTrackMov.superResolve(val2Use);
%% plot
frames = 1:100;

%MPTrackMov.showCorrLoc();

%% showFrame

%MPTrackMov.showFrame(80);
%MPTrackMov.showParticle;

%% tracking
trackParam.euDistXY = 600;
trackParam.euDistZ  = 600;
MPTrackMov.trackParticle(trackParam);
traces = MPTrackMov.getTraces;
%% plot
%MPTrackMov.showTraces;

%% eval accuracy
MPTrackMov.evalAccuracy

%% get RMSD

[RMSD,D] = MPTrackMov.getRMSD('2D');

%% Intensity and SNR
[int,SNR] = MPTrackMov.getAvgIntensity;

%% getTraces 3D
frame = 1:200;
idx2Trace = 5;
ROIradius = 20;
frameRate = 3;
scaleBar  = 500; %in nm 
%MPTrackMov.getTracesMovie
%MPTrackMov.getTraces3DMovie(frame,idx2Trace,frameRate);
MPTrackMov.getPartMovie(frame,idx2Trace,ROIradius,frameRate);
