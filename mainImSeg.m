clear 
close all
clc

%% User Input
prompt = {'Enter initial frame','Enter number of frames to analyze ',...
    'Enter the desired threshold sensitivity ([0 1])'};
dlgTitle = 'User input for image segmentation';
numLines = 1;
defaultVal = {'1','11','0.6'};
answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);

%% Checking user input
assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')

iniFra = round(str2double(answer(1)));
assert(~isnan(iniFra),'Initial frame should be numerical');%If not a number
nFrame = round(str2double(answer(2)));
assert(~isnan(nFrame),'Number of Frame should be numerical');%If not a number

finFra = iniFra + nFrame - 1;
assert(finFra>=iniFra, 'WTF not expected')

assert(iniFra>0, 'indexing starts at 1')

Threshold = str2double(answer(3));
assert(~isnan(Threshold),'Number of Frame should be numerical');%If not a number
assert(and(Threshold > 0, Threshold <=1),...
    'Threshold sensitivity should be between 0 and 1');

connectivity = 216; %3D connectivity for binarization
fileExt = '.tif'; %Extension of the files to be analyze
outputName  = 'SegmentedStacks';
%% Loading

[filename,path] = uigetfile('*.tif','Select a file to segment');
file2Load = [path filesep filename];
fileInfo    = Load.Movie.tif.getinfo(file2Load);
mkdir(path,outputName);
outDir = [path filesep outputName];
warning('off','all');
%Check number of Frame
tNframes = fileInfo.Frame_n;

if tNframes>finFra
    frames2load = iniFra:tNframes;
    warning('Requested number of frame is larger than the number of frame in the file, adapted the frame limit and run')

else
    frames2load = iniFra:finFra;
end

IM     = Load.Movie.tif.getframes(file2Load, frames2load); %Loading on of the frame
warning('on','all');

disp('DONE with loading --------------')

%% Processing
disp('Now doing 3D gauss filtering this can take about 3 minutes')
% size of gauss filter
S = 1;
% size of pixel in z vs x/y
pixZ  = 4;
zFactor = 2;
sigma = [S,S,S*zFactor/pixZ];
IMs = imgaussfilt3(IM, sigma);
disp('DONE with filtering ------------')

 diskDim = 4;
  [gBW,aBW] = imSegmentation.segmentStack(IMs,'threshold',Threshold,...
            'connectivity',connectivity,'diskDim',diskDim);

BWadapt  = aBW;
BWglobal = gBW;
disp('Storing global results')

tifName = [outDir filesep 'Seg_global_PolDark' filename];
dataStorage.BinaryTiff(tifName,BWglobal);
BWglobal = imcomplement(BWglobal);
tifName = [outDir filesep 'Seg_global_poreDark' filename];
dataStorage.BinaryTiff(tifName,BWglobal);

disp('Storing adaptive results')
tifName = [outDir filesep 'Seg_adapt_PolDark' filename];
dataStorage.BinaryTiff(tifName,BWadapt);
BWadapt = imcomplement(BWadapt);
tifName = [outDir filesep 'Seg_adapt_poreDark' filename];
dataStorage.BinaryTiff(tifName,BWadapt);

%Save useful information about how the binarization was done
infoFileName = [outDir filesep 'info.txt'];
fid = fopen(infoFileName,'wt');
fprintf(fid,'This file contains information intended for the user of segmentStack.m\n');
fprintf(fid,' In such a way that the user knows what variable value were used.\n\n');
fprintf(fid,'Adaptive Threshold sensitivity: %0.1f\n',Threshold);
fprintf(fid,'Number of frame analyzed: %d/%d\n',nFrame,tNframes);
fprintf(fid,'3D Gaussian filtering: S = %d; pixZ  = %d; sigma = [S,S,S*%d/pixZ]\n',S,pixZ,zFactor);
fprintf(fid,'Three-dimensional connectivity - BWareaopen: %d\n',connectivity);
fprintf(fid,'strel disk dimension for imopen: %d\n',diskDim);
fclose(fid);

%%
h = msgbox('The Data were succesfully saved !', 'Success');