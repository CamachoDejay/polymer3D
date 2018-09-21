%%
clc 
close all
clear 

path2ZCal = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - maxObjCorr';
path2SRCal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - maxObjCorrPSFE';

path2File = 'E:\Data\Leuven Data\2018\06-June\27\trackingCal - maxObjCorr\TrackingX\FluoBeads200_TrackingX_320ms__1';
path2Cal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - maxObjCorrPSFE\zStackFluoBeads200_S3_270618__1';

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

MPTrackMov.showCorrLoc();

%% showFrame

MPTrackMov.showFrame(80);
%MPTrackMov.showParticle;

%% tracking
trackParam.euDistXY = 300;
trackParam.euDistZ  = 300;
MPTrackMov.trackParticle(trackParam);
traces = MPTrackMov.getTraces;
%% plot
MPTrackMov.showTraces;

%% eval accuracy
MPTrackMov.evalAccuracy

%% Intensity and SNR
[int,SNR] = MPTrackMov.getAvgIntensity;

%% getTraces 3D
frame = 1:100;
idx2Trace = 5;
ROIradius = 10;
frameRate = 3;

%MPTrackMov.getTracesMovie
MPTrackMov.getTraces3DMovie(frame,idx2Trace,frameRate);
%MPTrackMov.getPartMovie
