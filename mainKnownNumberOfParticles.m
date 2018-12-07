clear 
clc
close all
%% User input
delta = 50;
nParticles = 2;
pxSize = 95;
minDist = 6; %in pixels
scaleBar = 2; %in um
tail = 20;
info.type = 'transmission';

toAnalze = 'folder';
outputFolder = 'Results';
%% Create the Movie object
switch toAnalze
    case '.mp4'
        [folder2Mov] = Core.Movie.getFileInPath(path2File, toAnalze);
        p2file = [folder2Mov.folder filesep folder2Mov.name];
                 
    case '.ome.tif'
        [folder2Mov,path,outDir] = Load.Folder(toAnalze,outputFolder);
        mkdir(outDir);
        
        
    case 'folder'
        [path] = uigetdir();
        tmp = dir(path);
        tmp = tmp(cell2mat({tmp.isdir}));
        tmp = tmp(3:end);
        folder2Mov = [];
        for i = 1:size(tmp,1)
            path2File = [tmp(i).folder filesep tmp(i).name];
            file2Analyze = Core.Movie.getFileInPath(path2File,'.ome.tif');
            
            folder2Mov = [folder2Mov file2Analyze];
        end
        outDir = [path filesep 'Results'];   
    otherwise
        error('Unknown extension');


end

assert(~isempty(folder2Mov),'Error no file was found, check that you put the correct analysis type');
%% detection of the center of the beads
allData = struct('fileName',[],'locPos',[]);
allData(size(folder2Mov,2)).locPos = [];

for i =1: size(folder2Mov,2)
    
    currentPath = folder2Mov(i).folder;
    switch toAnalyze
        case '.mp4'
            v = VideoReader(p2file);
            nFrames = floor(v.Duration*v.FrameRate);
            fullStackIn = zeros(v.Height,v.Width,nFrames);
            for j = 1:nFrames
                frame = readFrame(v);

                fullStackIn(:,:,j) = rgb2gray(frame);
            end
        otherwise
            myMov = Core.Movie(currentPath,info,'.tif');
            fullStack = myMov.getFrame(1);
            frame = fullStack.Cam1;
            if size(frame,2) > 400
                myMov.cropIm;
            end
            fullStack = myMov.getFrame;

            fullStackIn = imcomplement(fullStack.Cam1);
    end
            frame = 5;

    [pos] = goldProj.nMaxDetection (fullStackIn(:,:,frame),nParticles,minDist);

    x0 = pos(:,2);
    y0 = pos(:,1);

    pos = [y0 x0];
    cropPos = round(mean(pos,1));

    figure
    imagesc(fullStackIn(:,:,frame))
    hold on
    plot(pos(:,2),pos(:,1),'r+')

    %% Cropping Movie

    fullStackIn = fullStackIn(cropPos(1)-delta:cropPos(1)+delta, cropPos(2)-delta:cropPos(2)+delta,:);
    %% Data manipulation
    nFrames = size(fullStackIn,3);
    x = 1:size(fullStackIn,2);
    y = 1:size(fullStackIn,1);
    [domX,domY] = meshgrid(x,y);
    dom(:,:,1) = domX;
    dom(:,:,2) = domY;

    data2Store = zeros(nFrames,2,nParticles);
    fitMov = zeros(size(fullStackIn));
    h = waitbar(0,'Fitting Data');

    for j = 1:nFrames
        currentFrame = double(fullStackIn(:,:,j));

    [pos] = goldProj.nMaxDetection (currentFrame,nParticles,minDist);

        x0 = pos(:,2);
        y0 = pos(:,1);

        [gPar,resnorm,res] = Localization.Gauss.MultipleFitting(currentFrame,x0,y0,dom,nParticles); 
        F = Localization.Gauss.MultipleGauss(gPar, dom,nParticles);

        if j>1

            newOrder = goldProj.simpleTracking(gPar(5:end),prevPos);
            %reshaping and sorting with new order
            gPos = reshape(gPar(5:end),[2,nParticles]);
            gPos = gPos(:,newOrder);
            prevPos = reshape(gPos,[1,nParticles*2]);
            gPos = reshape(gPos,[1,2,nParticles]);
            data2Store(j,:,:) = gPos;
            clear gPos

        else

            gPos = reshape(gPar(5:end),[2,nParticles]);
            gPos = reshape(gPos,[1,2,nParticles]);
            data2Store(j,:,:) = gPos;
            %store the gPar in the prev for checking particle order
            prevPos = gPar(5:end);
            clear gPos
        end

        fitMov(:,:,j) = F;
        waitbar(j/nFrames,h,'Fitting Data')
    end
        
    filename = [currentPath filesep 'LocalizationData.mat'];
    save(filename,'data2Store');
    allData(i).locPos = data2Store;
    allData(i).fileName = currentPath;
    close(h);
    
    %% display figure
    data2plot = data2Store *pxSize;
    cm = [mean(mean(data2plot(:,1,:))) mean(mean(data2plot(:,2,:)))];
    Fig = figure;
    hold on
    for j = 1 : nParticles
        scatter(data2plot(:,1,j)-cm(1),data2plot(:,2,j)-cm(2),20,'filled')

    end
    axis image
    xlabel('X Position (px, 1px = 95nm)')
    ylabel('Y Position (px, 1px = 95nm)')
    title('All localized spot');
    filename = [currentPath filesep 'LocalizationDensity.fig'];
    saveas(Fig,filename);
    
    % Cropping full stack
switch toAnalze
    case '.mp4'
        fullStackNorm = fullStackIn;
    otherwise
        fullStackNorm = imcomplement(fullStackIn);
end
%% MovieMaker

frameRate = 30;
filename = [currentPath filesep 'TrackMovie.gif'];

goldProj.makeTraceMovie(data2Store,fullStackNorm,filename,frameRate,scaleBar,tail);

end

filename = [outDir filesep 'allLoc'];
save(filename,'allData');

%% Calculate correlation
% 
% %xcorr??
% euclDistX = sqrt((data2Store(:,1,1)- data2Store(:,1,2)).^2);
% avgMotX1 = data2Store(:,1,1) - mean(data2Store(:,1,1));
% avgMotX2 = data2Store(:,1,2) - mean(data2Store(:,1,2));
% avgMotX3 = data2Store(:,1,3) - mean(data2Store(:,1,3));
% avgMotX4 = data2Store(:,1,4) - mean(data2Store(:,1,4));
% 
% 
% %RXY = corrcoef(x,y);
% 
% %process MP4;
