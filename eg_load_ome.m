% the purpose of this program is to have a minimal working example of how
% to load a ".ome" file (micromanager) into matlab. We should be able to
% handle single and multi cam data, and extract all available information
% form the header.

clear
close all
clc

% at the moment this code is not great as I dont request the user for input
% regarding the file to load.

folder_path = 'path2folder';
main_file = 'name_of_main_file.ome.tif';
% note that an ome is just a moditiaction of a tif file.

path2data = [folder_path filesep main_file];

% getting information from header
[frameInfo, movieInfo, ~] = MovieLoad.ome.getomeinfo(path2data);
nFrames = size(frameInfo,2);

frames2load = 1:2;
[cam1, cam2, ~] = MovieLoad.ome.loadMulticamOME(frameInfo, movieInfo, frames2load);
