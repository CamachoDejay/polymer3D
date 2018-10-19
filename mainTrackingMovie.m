%%
clc 
close all
clear 

path2ZCal = 'M:\Data\Leuven Data\2018\06-June\27\ZCal - maxObjCorr';
path2SRCal = 'M:\Data\Leuven Data\2018\06-June\27\2DCal - maxObjCorrPSFE';

path2File = 'M:\Data\Leuven Data\2018\06-June\29\1K - 0_25mgmL\TL-FluoBeads200nm-PIC0_25mgmL-1K_2';
path2Cal = 'M:\Data\Leuven Data\2018\06-June\27\2DCal - maxObjCorrPSFE\zStackFluoBeads200_S3_270618__1';

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
trackParam.euDistXY = 1000;
trackParam.euDistZ  = 1300;
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
frame = 1:500;
idx2Trace = 1;
ROIradius = 15;
frameRate = 5;

%MPTrackMov.getTracesMovie(frame,idx2Trace,ROIradius,frameRate);
MPTrackMov.getTraces3DMovie(frame,idx2Trace,ROIradius,frameRate);
%MPTrackMov.getPartMovie
