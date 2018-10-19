%%
clc 
close all
clear 

path2ZCal = 'E:\Data\Leuven Data\2018\ZHao\181015 - Calibration\ZCal\Good Cal';
path2SRCal = [];

path2File = 'E:\Data\Leuven Data\2018\ZHao\181016 - Dry objective trapping\200nmGoldBeadsTrappingFoucsingLayer1_Initial_7';
path2Cal = 'E:\Data\Leuven Data\2018\ZHao\181015 - Calibration\2DCal\200nmFluoBeadsCalPSFE_3';
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