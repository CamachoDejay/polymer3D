% The aim of this code is to analyze a folder containing Z-stack with 
%structure e.g. porous material. 

% The program segments the stacks (e.g. binarization pore-material) and
% save the information for further analysis

clear
close all
clc

%% User Input
prompt = {'Enter number of frame to analyze ',...
    'Enter the desired threshold sensitivity ([0 1])'};
dlgTitle = 'User input for image segmentation';
numLines = 1;
defaultVal = {'244','0.6'};
answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);

%% Checking user input
assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')

nFrame = str2double(answer(1));
assert(~isnan(nFrame),'Number of Frame should be numerical');%If not a number

Threshold = str2double(answer(2));
assert(~isnan(Threshold),'Number of Frame should be numerical');%If not a number
assert(and(Threshold > 0, Threshold <=1),...
    'Threshold sensitivity should be between 0 and 1');

connectivity = 216; %3D connectivity for binarization
fileExt = '.tif'; %Extension of the files to be analyze
outputName  = 'SegmentedStacks';
%% Loading Data
%Load folder, and create a folder for data output.
[file2Analyze,currentFolderName,outDir] = Load.Folder(fileExt,outputName);
assert(~isempty(file2Analyze), sprintf('no %s found in the directory', fileExt));
%% load full stack
nFiles = size(file2Analyze,1);

for i = 1 : nFiles

    disp(['Loading stack --------------' file2Analyze(i).name])
    disp('This can take a few minutes ~2')
    tic
    path2Stacks = strcat(file2Analyze(i).folder,filesep);
    tmpName = file2Analyze(i).name;
    p2file      = strcat(path2Stacks,tmpName);
    fileInfo    = Load.Movie.tif.getinfo(p2file);
    warning('off','all');
    %Check number of Frame
    tNframes = fileInfo.Frame_n;
    assert(tNframes>=nFrame,'Requested number of frame is larger than the number of frame in the file')
    
    IM     = Load.Movie.tif.getframes(p2file, 1:nFrame); %Loading on of the frame
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
    zFactor = 3;
    sigma = [S,S,S*zFactor/pixZ];
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
    h = waitbar(0, sprintf('Starting segmentation of stack %d/%d...',i,size(file2Analyze,2)));
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
        waitbar(j/lastIDX,h,sprintf('Segmentation of stack %d/%d... %d percent achieved',...
            i,size(file2Analyze,2), round(j/lastIDX*100)));
    end

    %%%%%%%%%%%%%%% Data storing %%%%%%%%%%%%%%%
    
    % now we save segmented images
    disp('Storing global results')
    tifName = [outDir filesep 'Seg_global_' tmpName];
    dataStorage.BinaryTiff(tifName,BWglobal);

    disp('Storing adaptive results')
    tifName = [outDir filesep 'Seg_adapt_' tmpName];
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
    fprintf(fid,'Image partition in chunck: %d x %d x %d \n',dIDX,dIDX,nFrame);
  
    
    fclose(fid);
    close(h);
end

%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%

figure(1)
shg
SE = strel('disk',3);

for j = round(linspace(1,nFrame,10))
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

