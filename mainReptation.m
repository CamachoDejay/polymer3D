%% main Reptation movie
clc;
clear;
clc;
%% User Input
% path to the callibration
file.path =  'D:\Dropbox\JohannesDataTMP\dna_50msec OD 1.9__1';
path2Cal  =  ['D:\Dropbox\JohannesDataTMP\2Dcali'];
file.ext = '.ome.tif';

%% Info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.frame2Load = 'all'; % 'all' or a range of number e.g. 1:100
info.fitMethod  = 'Phasor'; %Phasor or Gauss (need to be the same as ZCal if using PSFE
info.zMethod = 'Intensity'; %Intensity, 3DFit or PSFE
info.calibrate = false; %true to recalibrate;

%% Data Loading

myMovie = Core.ReptationMovie(file,path2Cal,info);

myMovie.calibrate;
myMovie.giveInfo;

%% get a molecule
myMovie.showFrame(1,5)


allAxesInFigure = findall(figure(1),'type','axes');

F = figure(2);
hNew = copyobj(allAxesInFigure(6),F);
set(hNew, 'pos', [0.23162 0.2233 0.72058 0.63107])
h = drawrectangle();

pos = h.Position;

data = myMovie.cropAllFrames(pos);

%% play around with data

