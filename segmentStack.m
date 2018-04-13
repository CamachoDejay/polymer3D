% program used to segment fulls stacks and save information for further
% analysis
clear
close all
clc

%% User Input
Threshold = 0.6; %number between 0 and 1 (% of max intensity)
% used during testing 2, normal should 244
nFrame = 244; %n of frame to analyze
fileExt = '.tif'; %Extension of the files to be analyze
outputName  = 'SegmentedStacks';
%% Loading Data
%Load folder, and create a folder for data output.
[file2Analyze,currentFolderName,outDir] = Misc.loadFolder(fileExt,outputName);

%% load full stack
%TODO: Loop through stack
stack2Load  = 1;

disp(['Loading stack --------------' file2Analyze(stack2Load).name])
disp('This can take a few minutes ~2')
tic
path2Stacks = strcat(file2Analyze(stack2Load).folder,filesep);
tmpName = file2Analyze(stack2Load).name;
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
    
    %% Store data
    % now we save segmented images
    disp('Storing global results')
    % Save global T
    tifName = [outDir filesep 'Seg_global_' tmpName];
    dataStorage.saveBinaryTiff(tifName,BWglobal);

    disp('Storing adaptive results')
    tifName = [outDir filesep 'Seg_adapt_' tmpName];
    dataStorage.saveBinaryTiff(tifName,BWglobal);

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


