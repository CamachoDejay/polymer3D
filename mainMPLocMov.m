%%
clc 
close all
clear 

path2ZCal = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\ZCal - NormObjCorr';
path2SRCal = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\2DCal - normObjCorrPSFE';

path2File = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\ZCal - NormObjCorr\zStackFluoBeads200_PIC1_270618__2';
path2Cal = 'D:\Documents\Unif\PhD\2018-Data\06 - Jun\27\ZCal - NormObjCorr\zStackFluoBeads200_PIC1_270618__1';

detectParam.delta = 6;
detectParam.chi2 = 80;
%%
calib = Core.MPCalibration(path2Cal);

%%

MPLocMov = Core.MPLocMovie(path2File,calib.getCal,path2SRCal,path2ZCal);

%% 
MPLocMov.showFrame(40);

%% Detection

MPLocMov.giveInfo
%find candidate
MPLocMov.findCandidatePos(detectParam);

%fit position
MPLocMov.SRLocalizeCandidate;

%% Data correction
rot = true;
refPlane = 5;
MPLocMov.applySRCal(rot,refPlane);
%% e-Z transformation
MPLocMov.applyZCal;

%% Plane consolidation

MPLocMov.consolidatePlanes
%% Super resolve

MPLocMov.superResolve;
%% plot

MPLocMov.showCorrLoc;