%%
clc 
close all
clear 

path2ZCal = 'D:\Documents\Unif\PhD\2019-Data\04 - Apr\extended\test\ZCal';
path2SRCal = 'D:\Documents\Unif\PhD\2019-Data\04 - Apr\extended\test\2DCal';

path2File = 'D:\Documents\Unif\PhD\2019-Data\04 - Apr\extended\test\XTrack';
path2Cal = 'D:\Documents\Unif\PhD\2019-Data\04 - Apr\extended\test\2DCal';
info.type = 'transmission';
info.runMethod = 'load';
detectParam.delta = 6;
detectParam.chi2 = 80;
%%
calib = Core.MPCalibration(path2Cal);

%%

MPLocMov = Core.MPLocMovie(path2File,calib.getCal,info,path2SRCal,path2ZCal);

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