%% 
%The aim of this code is to receive a z-stack/ a serie of z-stack and to
%calculate the relation between the z-position and the ellipticity to serva
%as a calibration in 3D particle tracking.

% It is currently written to take .mat file (the file are prealably cropped
% as we do no consider multiplane data yet.

clear;
clc;
close all;

%% User input 
FWHM_nm = 350;%in nm
pxSize  = 105;%in nm
szWindow = 6;
zSpacing = 50; %in nm

setupInfo.FWHM = FWHM_nm;
setupInfo.pxSize = pxSize;
setupInfo.szWindow = szWindow;
setupInfo.zSpacing = zSpacing;
%% Loading of the data
[fileName,folder]=uigetfile({'*.mat'},'Select a file to Open');
[~,fileTif,~]=fileparts(fileName);
[~,fileOME,~]=fileparts(fileTif);

path2data = [folder fileName];

tmp=load(path2data);
name=fields(tmp);
imStack=tmp.(name{1});

%% Z-Calibration
[zAxis,ellipAxis] = zCalibration(setupInfo,imStack);

figure
scatter(zAxis,ellipAxis);



