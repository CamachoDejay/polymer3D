%%%%%%%%%%%%%%%%%%%%%% PLANE CALIBRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: Folder containing folder containing .ome.tif file that are 2D    %
% calibration.                                                            %
%                                                                         %
% Results: 2D Calibration will be saved (ROI + position of planes in z    %
%                                                                         %
% if the offTarget returned to the command window is >20 nm there might be%
% something wrong with either the analysis,the dataset or the alignment of%
% the camera                                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

file.path  = 'D:\Documents\Unif\PhD\2020-Data\08 - August\Johannes\Movies boris\2Dcal';

file.ext   = '.ome.tif';
info.runMethod = 'load';
info.nChan = 4;
%% 
calib = Core.MPPlaneCalibration(file,info);

calib.retrieveMovies;
calib.calcIndivCal; 
calib.calcCombinedCal;

calib.showCal(1)
calib.offTarget;
calib.save;
