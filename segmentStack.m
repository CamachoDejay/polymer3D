% The aim of this code is to analyze a folder containing Z-stack with 
%structure e.g. porous material. 

% The program segments the stacks (e.g. binarization pore-material) and
% save the information for further analysis

clear
close all
clc

%% User Input
Threshold = 0.6; %number between 0 and 1 (% of max intensity)
% used during testing 2, normal should 244
nFrame = 244; %n of frame to analyze
connectivity = 216; %3D connectivity for binarization

fileExt = '.tif'; %Extension of the files to be analyze
outputName  = 'SegmentedStacks';
%% Loading Data
%Load folder, and create a folder for data output.
[file2Analyze,currentFolderName,outDir] = Misc.loadFolder(fileExt,outputName);

%% load full stack

for i = 1 : size(file2Analyze,2)

    disp(['Loading stack --------------' file2Analyze(i).name])
    disp('This can take a few minutes ~2')
    tic
    path2Stacks = strcat(file2Analyze(i).folder,filesep);
    tmpName = file2Analyze(i).name;
    p2file      = strcat(path2Stacks,tmpName);
    fileInfo    = loadMovie.tif.getinfo(p2file);
    warning('off','all');
    IM     = loadMovie.tif.getframes(p2file, 1:nFrame); %Loading on of the frame
    warning('on','all');
    toc
    disp('DONE with loading --------------')

    %%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%
    
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

    %%%%%%%%%%%%%%% Segmenting %%%%%%%%%%%%%%%
    
    disp('Now doing segmentation this can take a few minutes ~10')
    imSize = size(IMs);
    BWglobal = false(imSize);
    BWadapt  = false(imSize);

    % get list of indices
    dIDX = 512; % to split stack in 512*512*Z chunks
    xi = 1:dIDX:imSize(1);
    xf = dIDX:dIDX:imSize(1);
    IDX = zeros(length(xi)*length(xf),4);

    [XX, YY] = meshgrid(xi,xf);
    IDX(:,1) = XX(:);
    IDX(:,4) = YY(:);
    [XX, YY] = meshgrid(xf,xi);
    IDX(:,2) = XX(:);
    IDX(:,3) = YY(:);
    clear XX YY xi xf 

    lastIDX = size(IDX,1);
    lastStr = num2str(lastIDX);
    diskDim = 4;
    for j = 1:lastIDX
        tic
        iStr = num2str(j);
        xi = IDX(j,1);
        xf = IDX(j,2);
        yi = IDX(j,3);
        yf = IDX(j,4);
        tmp = IMs(xi:xf,yi:yf,:);

        th = adaptthresh(tmp,Threshold,'neigh',[301 301 151],'Fore','bright');
        BW = imbinarize(tmp,th);
        BW = ~BW;
        BW = bwareaopen(BW,connectivity);
        SE = strel('disk',diskDim);
        BW = imopen(BW,SE);
        BWadapt (xi:xf,yi:yf,:) = BW;
        disp(['Done adaptive step ' iStr '/' lastStr ])

        BW = imbinarize(tmp,'global');
        BW = ~BW;
        BW = bwareaopen(BW,connectivity);
        SE = strel('disk',diskDim);
        BW = imopen(BW,SE);
        BWglobal(xi:xf,yi:yf,:) = BW;
        disp(['Done global step ' iStr '/' lastStr ])
        toc
    end

    %%%%%%%%%%%%%%% Data storing %%%%%%%%%%%%%%%
    
    % now we save segmented images
    disp('Storing global results')
    tifName = [outDir filesep 'Seg_global_' tmpName];
    dataStorage.saveBinaryTiff(tifName,BWglobal);

    disp('Storing adaptive results')
    tifName = [outDir filesep 'Seg_adapt_' tmpName];
    dataStorage.saveBinaryTiff(tifName,BWglobal);
    
    %Save useful information about how the binarization was done
    infoFileName = [outDir filesep 'info.txt'];
    fid = fopen(infoFileName,'wt');
    fprintf(fid,'This file contains information intended for the user of segmentStack.m\n');
    fprintf(fid,' In such a way that the user knows what variable value were used.\n\n');
    fprintf(fid,'Adaptive Threshold sensitivity: %0.1f\n',Threshold);
    fprintf(fid,'Three-dimensional connectivity - BWareaopen: %d\n',connectivity);
    fprintf(fid,'strel disk dimension for imopen: %d\n',diskDim);
    fprintf(fid,'Image partition in chunck: %d x %d x %d \n',dIDX,dIDX,nFrame);
    
    fclose(fid);

    %%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%
    if i==1 %only plot for the first stack
        figure(1)
        shg
        SE = strel('disk',3);
        for j = 1:50
            subplot(1,2,1)
            A = IMs(:,:,j);
            B = BWglobal(:,:,j);
            B = bwperim(B);
            B = imdilate(B,SE);
            C = imfuse(A,B,'ColorChannels',[2 1 0]);
            imagesc(C)
            axis image
            title('global')

            subplot(1,2,2)
            A = IMs(:,:,j);
            B = BWadapt(:,:,j);
            B = bwperim(B);
            B = imdilate(B,SE);
            C = imfuse(A,B,'ColorChannels',[2 1 0]);
            imagesc(C)
            axis image
            title('adaptive')
            waitforbuttonpress

        end
    end
end

