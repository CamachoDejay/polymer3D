% the purpose of this program is to have a minimal working example of how
% to load a ".ome" file (micromanager) into matlab. We should be able to
% handle single and multi cam data, and extract all available information
% form the header.

clear
close all
clc

% at the moment this code is not great as I dont request the user for input
% regarding the file to load.
[fileName,folder]=uigetfile({'*.tif'},'Select a file to Open');
[~,fileTif,~]=fileparts(fileName);
[~,fileOME,~]=fileparts(fileTif);

path2data = [folder fileName];

% getting information from header
[frameInfo, movieInfo, ~] = loadMovie.ome.getInfo(path2data);


frames2load = 1:min(movieInfo.maxFrame);%Use min in case the number of frame
%on the two camera is not equal.
[cam1, cam2, ~] = loadMovie.ome.load(frameInfo, movieInfo, frames2load);
