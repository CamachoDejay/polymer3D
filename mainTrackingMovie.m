%%
clc 
close all
clear 

path2ZCal = 'E:\Data\Leuven Data\2018\ZHao\181015 - Calibration\ZCal\Good Cal';
path2SRCal = [];

path2File = 'E:\Data\Leuven Data\2018\ZHao\181018 - 400nm GoldBeads Trapping\GoldBeads400nm - Water\GoldBeads400nmTransmission_IntermediateStage_1';
path2Cal = 'E:\Data\Leuven Data\2018\ZHao\181015 - Calibration\2DCal\200nmFluoBeadsCalPSFE_3';

detectParam.delta = 6;
detectParam.chi2 = 40;
info.runMethod = 'run';
info.type = 'transmission';
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
frames = 1:500;

MPTrackMov.showCorrLoc(frames);

%% showFrame

%MPTrackMov.showFrame(80,5);
%MPTrackMov.showParticle;

%% tracking
trackParam.euDistXY = 400;
trackParam.euDistZ  = 400;
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
