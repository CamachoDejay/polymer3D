function checkSeg (path, idx)
% Extract Tif file from the given path and do a couple of assertion before
% loading anything
folder2Segmentation = [path filesep 'SegmentedStacks'];
assert(ischar(path),'Path given should be a char');
assert(isfolder(path),'The path given is not a folder');
assert(isfolder(folder2Segmentation),'No segmented stack found in the given directory');
fileExt = '.tif';
%Extract the part of the folder that is a tif file
folderContent = dir(path);
index2Images   = contains({folderContent.name},fileExt);
file2Analyze = folderContent(index2Images);
assert(~isempty(file2Analyze),'No file found in the selected directory');
assert(idx<=length(file2Analyze),'Index requested exceeds number of file in the directory');

%Here we load the raw data and do the filtering as done in the segmentation
%procedure
disp('Loading the raw data')
disp(['Loading stack --------------' file2Analyze(idx).name])
path2Stacks = strcat(file2Analyze(idx).folder,filesep);
tmpName = file2Analyze(idx).name;
p2file      = strcat(path2Stacks,tmpName);
fileInfo    = Load.Movie.tif.getinfo(p2file);
warning('off','all');


frames = fileInfo.Frame_n;

IM     = Load.Movie.tif.getframes(p2file, 1:frames);

disp('DONE with loading --------------')

%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%
% IMtest = IM(1:500,1:500,1:100);
disp('Now doing 3D gauss filtering this can take about 3 minutes')

tic
% size of gauss filter
S = 1;
% size of pixel in z vs x/y
pixZ  = 5;
zFactor = 3;
sigma = [S,S,S*zFactor/pixZ];
IMs = imgaussfilt3(IM, sigma);
toc
disp('DONE with filtering ------------')

%Here we load the segmented stack that correspond to the selected file
disp('Loading the segmented data')

folderContent = dir(folder2Segmentation);
index2Images   = contains({folderContent.name},fileExt);
fileAnalyze = folderContent(index2Images);
assert(~isempty(fileAnalyze),'No file found in the selected directory');

index2Stack = contains({fileAnalyze.name},file2Analyze(idx).name);
correspSegmentedStack = fileAnalyze(index2Stack);
assert(~isempty(correspSegmentedStack),'No corresponding segmented stack found in the directory');
assert(length(correspSegmentedStack)==2,'One of the segmented stack is missing (global or adaptive) please check');

path2Adapt = strcat(correspSegmentedStack(1).folder,filesep);
tmpName = correspSegmentedStack(1).name;
p2file  = strcat(path2Adapt,tmpName);
fileInfo= Load.Movie.tif.getinfo(p2file);
warning('off','all');

frames = fileInfo.Frame_n;

BWadapt     = Load.Movie.tif.getframes(p2file, 1:frames);

path2Global = strcat(correspSegmentedStack(2).folder,filesep);
tmpName = correspSegmentedStack(2).name;
p2file  = strcat(path2Global,tmpName);
fileInfo= Load.Movie.tif.getinfo(p2file);
warning('off','all');

frames = fileInfo.Frame_n;

BWglobal     = Load.Movie.tif.getframes(p2file, 1:frames);

%Plot to give an idea about the goodness of the segmentation compare to the
%raw data

figure(1)
shg
SE = strel('disk',3);

for j = round(linspace(1,frames,10))
    subplot(1,2,1)
    A = IMs(:,:,j);
    B = BWglobal(:,:,j);
    B = bwperim(B);
    B = imdilate(B,SE);
    C = imfuse(A,B,'ColorChannels',[2 1 0]);
    imagesc(C)
    axis image
    title(['global - frame ' num2str(j)] )

    subplot(1,2,2)
    A = IMs(:,:,j);
    B = BWadapt(:,:,j);
    B = bwperim(B);
    B = imdilate(B,SE);
    C = imfuse(A,B,'ColorChannels',[2 1 0]);
    imagesc(C)
    axis image
    title(['adaptive - frame ' num2str(j)])
    waitforbuttonpress

end