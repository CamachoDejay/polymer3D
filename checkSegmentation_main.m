clear
close all
clc

idx2File = 1;
pathForTest = 'N:\DATA\Leica\PhD\2018\10062018\tif cold 1';
path = uigetdir;

%%
checkSegmentation(path, idx2File);