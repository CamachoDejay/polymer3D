%%
clc 
close all
clear 

path2ZCal = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - maxObjCorr';
path2SRCal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - maxObjCorrPSFE';

path2File = 'E:\Data\Leuven Data\2018\06-June\29\1K - 0_25mgmL\TL-FluoBeads200nm-PIC0_25mgmL-1K_2';
path2Cal = 'E:\Data\Leuven Data\2018\06-June\29\2DCal\zStackFLuoBeads_2D3DS3__1';

detectParam.delta = 6;
detectParam.chi2 = 80;
info.runMethod = 'load';
info.type = 'normal';
%%
calib = Core.MPCalibration(path2Cal);

%%

MPTrackMov = Core.MPTrackingMovie(path2File,calib.getCal,info,path2SRCal,path2ZCal);

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
frames = 1:100;

%MPTrackMov.showCorrLoc();

%% showFrame

%MPTrackMov.showFrame(80);
%MPTrackMov.showParticle;

%% tracking
trackParam.euDistXY = 1000;
trackParam.euDistZ  = 1000;
MPTrackMov.trackParticle(trackParam);
traces = MPTrackMov.getTraces;
%% plot
%MPTrackMov.showTraces;

%% eval accuracy
MPTrackMov.evalAccuracy

%% get RMSD

[RMSD,D] = MPTrackMov.getRMSD('3D');

%% Intensity and SNR
[int,SNR] = MPTrackMov.getAvgIntensity;

%% getTraces 3D
frame =1:100;
idx2Trace = 1;
ROIradius = 12;
frameRate = 5;
scaleBar  = 500; %in nm 
MPTrackMov.getTracesMovie(frame,idx2Trace,ROIradius,frameRate,scaleBar);
%MPTrackMov.getTraces3DMovie(frame,idx2Trace,ROIradius,frameRate);
%MPTrackMov.getPartMovie(frame,idx2Trace,ROIradius,frameRate);
