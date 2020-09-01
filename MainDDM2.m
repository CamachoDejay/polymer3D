file.ext   = '.ome.tif';
path2Cal =  ['F:\3-9-2020Cal2D'];
dimension = '3D';
file.path  = ['F:\DataDDMPolystyrene\500nm\500nmI1in5Dil-BrightField-DDM_10ms\500nm1in5_3' ];

%% Storing info about the file
info.type = 'normal'; %normal or transmission
info.runMethod = 'load'; % load or run
info.calibrate = false; %true to recalibrate;


%% Input For DDM 
DDMInfo.nFrames =200;% number of frames to load into memory
DDMInfo.FramesToAnalyze = 120; %Number of frames to analyze
DDMInfo.PixelSize = 0.095;

DDMMovie = Core.DDMMovie(file,path2Cal,info,DDMInfo);
DDMMovie.calibrate();
DDMMovie.LoadAllFrames();
DDMMovie.main();

%% Fitting DDM Result
ToFit = Core.FitDDMMovie(DDMMovie.DDMOutput,'a*(1-exp(-x/b))+c');          
ToFit.FitDDMData;
ToFit.ShowFitResult;