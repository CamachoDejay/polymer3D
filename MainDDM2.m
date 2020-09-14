file.ext   = '.tif';
path2Cal =  [];
dimension = '3D';
file.path  = ['D:\Documents\Unif\PhD\2020-Data\09 - Sep\simBlinking' ];
correctDrift = false;


%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.calibrate = false; %true to recalibrate;

%% Input For DDM 
DDMInfo.nFrames = 500;% number of frames to load into memory
DDMInfo.FramesToAnalyze = 200 ; %Number of frames to analyze
DDMInfo.PixelSize = 0.2;

DDMMovie = Core.DDMMovie(file,path2Cal,info,DDMInfo);
DDMMovie.calibrate();
DDMMovie.LoadAllFrames(correctDrift);

%% Extract DDM signal 
tic 
profile('on')
DDMMovie.main('NumBins',50);
toc
profile('viewer')
%% Fitting DDM Result
ToFit = Core.FitDDMMovie(DDMMovie.DDMOutput,'a*(1-exp(-x/b))+c');          
ToFit.FitDDMData;
ToFit.ShowFitResult;