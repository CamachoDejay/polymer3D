% program used to segment fulls stacks and save information for further
% analysis

clear
close all
clc

%% User Input
Threshold = 0.6; %number between 0 and 1 (% of max intensity)
% used during testing 2, normal should 244
nFrame = 244; %n of frame to analyze
%% Loading Data
mainFolderName = uigetdir;
assert(ischar(mainFolderName),'User canceled the selection of file, excecution aborted');
%extract the name of the current folder
idx = strfind(mainFolderName,filesep) ;
currentFolderName = mainFolderName(idx(end)+1:end) ;
%Remove dots from the name
currentFolderName = regexprep(currentFolderName,'\.','_');

% generate folder to store output
segDir = [mainFolderName filesep 'SegementedStacks'];
status = mkdir(segDir);

%Extract the part of the folder that is a tif file
Folder_Content = dir(mainFolderName);
index2Images   = contains({Folder_Content.name},'.tif');
images2Analyze = Folder_Content(index2Images);
%% load full stack

stack2Load  = 1;

disp(['Loading stack --------------' images2Analyze(stack2Load).name])
disp('This can take a few minutes ~2')
tic
path2Stacks = strcat(images2Analyze(stack2Load).folder,filesep);
tmpName = images2Analyze(stack2Load).name;
p2file      = strcat(path2Stacks,tmpName);
fileInfo    = loadMovie.tif.getinfo(p2file);
warning('off','all');
IM     = loadMovie.tif.getframes(p2file, 1:nFrame); %Loading on of the frame
warning('on','all');
toc
disp('DONE with loading --------------')

%% Filtering
% IMtest = IM(1:500,1:500,1:100);
disp('Now doing 3D gauss filtering this can take about 3 minutes')
tic

% size of gauss filter
S = 1;
% size of pixel in z vs x/y
pixZ  = 5;
sigma = [S,S,S*3/pixZ];
IMs = imgaussfilt3(IM, sigma);
toc
disp('DONE with filtering ------------')

% for i = 1:nFrame
%     figure(1)
%     subplot(1,2,1)
%     imagesc(IM(:,:,i))
%     axis image
% 
%     subplot(1,2,2)
%     imagesc(IMs(:,:,i))
%     axis image
%     
%     waitforbuttonpress
% end
%% Segmenting
disp('Now doing segmentation this can take a few minutes ~10')
imSize = size(IMs);
BWglobal = false(imSize);
BWadapt  = false(imSize);

%test on a single frame - does not really correlate to what we see after
%usignt he 3D
% % % % imTest = IMs(:,:,10);
% % % % 
% % % % BW = imbinarize(imTest,'global');
% % % % BW = ~BW;
% % % % BW = bwareaopen(BW,216);
% % % % SE = strel('disk',4);
% % % % BWgTest = imopen(BW,SE);
% % % % 
% % % % th = adaptthresh(imTest,Threshold,'neigh',[251 251],'Fore','bright');
% % % % BW = imbinarize(imTest,th);
% % % % BW = ~BW;
% % % % BW = bwareaopen(BW,216);
% % % % SE = strel('disk',4);
% % % % BWaTest = imopen(BW,SE);
% % % % 
% % % % figure(1)
% % % % shg
% % % % SE = strel('disk',3);
% % % % subplot(1,2,1)
% % % % A = imTest;
% % % % B = BWgTest;
% % % % B = bwperim(B);
% % % % B = imdilate(B,SE);
% % % % C = imfuse(A,B,'ColorChannels',[2 1 0]);
% % % % imagesc(C)
% % % % axis image
% % % % title('global')
% % % % 
% % % % subplot(1,2,2)
% % % % A = imTest;
% % % % B = BWaTest;
% % % % B = bwperim(B);
% % % % B = imdilate(B,SE);
% % % % C = imfuse(A,B,'ColorChannels',[2 1 0]);
% % % % imagesc(C)
% % % % axis image
% % % % title('adaptive')

% get list of indices
dIDX = 512;
xi = 1:dIDX:imSize(1);
xf = dIDX:dIDX:imSize(1);
IDX = zeros(length(xi)*length(xf),4);

[XX, YY] = meshgrid(xi,xf);
IDX(:,1) = XX(:);
IDX(:,4) = YY(:);
[XX, YY] = meshgrid(xf,xi);
IDX(:,2) = XX(:);
IDX(:,3) = YY(:);
clear XX YY xi xf dIDX

lastIDX = size(IDX,1);
lastStr = num2str(lastIDX);

for i = 1:lastIDX
    tic
    iStr = num2str(i);
    xi = IDX(i,1);
    xf = IDX(i,2);
    yi = IDX(i,3);
    yf = IDX(i,4);
    tmp = IMs(xi:xf,yi:yf,:);
    
    th = adaptthresh(tmp,Threshold,'neigh',[301 301 151],'Fore','bright');
    BW = imbinarize(tmp,th);
    BW = ~BW;
    BW = bwareaopen(BW,216);
    SE = strel('disk',4);
    BW = imopen(BW,SE);
    BWadapt (xi:xf,yi:yf,:) = BW;
    disp(['Done adaptive step ' iStr '/' lastStr ])
    
    BW = imbinarize(tmp,'global');
    BW = ~BW;
    BW = bwareaopen(BW,216);
    SE = strel('disk',4);
    BW = imopen(BW,SE);
    BWglobal(xi:xf,yi:yf,:) = BW;
    disp(['Done global step ' iStr '/' lastStr ])
    toc
    
end


% figure(2)
% for i = 1:100;
% subplot(1,2,1)
% imagesc(IMs(1:512,1:2*512,i)); axis image
% 
% subplot(1,2,2)
% imagesc(BWadapt(1:512,1:2*512,i)); axis image
% 
% waitforbuttonpress
% end
%% Store data

disp('Storing global results')
% now we save segmented images
tifName = [segDir filesep 'Seg_global_' tmpName];
t = Tiff(tifName, 'w');
setTag(t,'ImageLength',size(BWglobal,1))
setTag(t,'ImageWidth',size(BWglobal,2))
setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
setTag(t,'BitsPerSample',1)
setTag(t,'SamplesPerPixel',1)
setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
t.write(BWglobal(:,:,1))

for i = 2:size(BWglobal,3)
%     disp(i)
    t.writeDirectory
    setTag(t,'ImageLength',size(BWglobal,1))
    setTag(t,'ImageWidth',size(BWglobal,2))
    setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
    setTag(t,'BitsPerSample',1)
    setTag(t,'SamplesPerPixel',1)
    setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
    t.write(BWglobal(:,:,i))
end
t.close 


disp('Storing adaptive results')
tifName = [segDir filesep 'Seg_adapt_' tmpName];
t = Tiff(tifName, 'w');
setTag(t,'ImageLength',size(BWadapt,1))
setTag(t,'ImageWidth',size(BWadapt,2))
setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
setTag(t,'BitsPerSample',1)
setTag(t,'SamplesPerPixel',1)
setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
t.write(BWadapt(:,:,1))

for i = 2:size(BWadapt,3)
%     disp(i)
    t.writeDirectory
    setTag(t,'ImageLength',size(BWadapt,1))
    setTag(t,'ImageWidth',size(BWadapt,2))
    setTag(t,'Photometric',Tiff.Photometric.MinIsBlack)
    setTag(t,'BitsPerSample',1)
    setTag(t,'SamplesPerPixel',1)
    setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky)
    t.write(BWadapt(:,:,i))
end
t.close
    
%%
figure(1)
shg
SE = strel('disk',3);
for i = 1:50
    subplot(1,2,1)
    A = IMs(:,:,i);
    B = BWglobal(:,:,i);
    B = bwperim(B);
    B = imdilate(B,SE);
    C = imfuse(A,B,'ColorChannels',[2 1 0]);
    imagesc(C)
    axis image
    title('global')
    
    subplot(1,2,2)
    A = IMs(:,:,i);
    B = BWadapt(:,:,i);
    B = bwperim(B);
    B = imdilate(B,SE);
    C = imfuse(A,B,'ColorChannels',[2 1 0]);
    imagesc(C)
    axis image
    title('adaptive')
    
    waitforbuttonpress
    
end


