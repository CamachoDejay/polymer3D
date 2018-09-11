%%
clc 
close all
clear 

path2ZCal = 'E:\Data\Leuven Data\2018\03-Mar\22\ZCalibration\';
path2SRCal = 'E:\Data\Leuven Data\2018\03-Mar\22\PlaneCaibration\';

path2File = 'E:\Data\Leuven Data\2018\03-Mar\22\ZCalibration\BeadsZCalibration_2';
path2Cal = '..\data\Multiplane\PlaneCalib\BeadsCalibrationZStack_1';

detectParam.delta = 6;
detectParam.chi2 = 80;
%%
calib = Core.MPCalibration(path2Cal);

%%

MPLocMov = Core.MPLocMovie(path2File,calib.getCal,path2SRCal,path2ZCal);

%% Detection

MPLocMov.giveInfo
%find candidate
MPLocMov.findCandidatePos(detectParam);

%fit position
MPLocMov.SRLocalizeCandidate;

%% Data correction
rot = false;
refPlane = 5;
MPLocMov.applySRCal(rot,refPlane);
%% e-Z transformation
MPLocMov.applyZCal;

%% Plane consolidation

MPLocMov.consolidatePlanes
%% STOPPED HERE NEED TO CODE SUPER RESOLVE POSITION SO WE GET 1 position per particle per time point !!! 

%% plot

MPLocMov.showCorrLoc;