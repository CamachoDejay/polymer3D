%%
clc;
close all;
clear;

path2ZCal = [];
path2SRCal = [];

path2File = 'F:\Data\Leuven Data\2019\05 - May\10\2D_cal_NOPSFE_1';
path2Cal  =  'E:\leuven data\05 May\13\2D Cal';


detectParam.delta = 6;
detectParam.chi2 = 20;
roiRad = 6;
info.runMethod = 'run';
info.type = 'normal';
info.zMethod = 'Intensity';
info.fitMethod = 'Phasor';
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
trackParam.euDistXY = 400;
trackParam.euDistZ  = 300;
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
