%%
clc 
close all
clear 

path2ZCal = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - NormObjCorr';
path2SRCal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE';

path2File = 'E:\Data\Leuven Data\2018\06-June\27\ZCal - NormObjCorr\zStackFluoBeads200_PIC1_270618__5';
path2Cal = 'E:\Data\Leuven Data\2018\06-June\27\2DCal - normObjCorrPSFE\zStackFluoBeads200_S3_270618__1';

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
trackParam.euDistZ  = 400;
MPTrackMov.trackParticle(trackParam);
%% plot
MPTrackMov.showTraces;

%% eval accuracy
MPTrackMov.evalAccuracy

%% Susana's figure


