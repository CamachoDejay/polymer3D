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

fr = data(:,:,:,1);
%fr = imgaussfilt3(fr,[2 2 1]);
fr = fr./median(median(fr,1),2);

BW = fr>1.2;

BW = bwareaopen(BW,100);

figure
imagesc(BW(:,:,4))

I = regionprops3(BW,'Volume','VoxelIDXList');

newImage = zeros(size(BW));

[val,idx] = max([I.Volume]);
pxList = I.VoxelIdxList{idx};

newImage(pxList) = 1;

figure
for i = 1:size(newImage,3)
   subplot(2,4,i)
   imagesc(newImage(:,:,i))
   axis image
end

%clean data
kernel = ones(3)/4;
for i = 1:size(newImage,3)
   cBW = newImage(:,:,8);
   cBW = bwareaopen(cBW,4);
   cBW = imfill(cBW,'holes');
   
   blurryImage = conv2(single(cBW),kernel, 'same');
   binary    = blurryImage >0.5;
    
end


