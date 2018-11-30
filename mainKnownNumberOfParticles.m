clear 
clc
close all
%% User input
delta = 30;
nParticles = 4;
pxSize = 95;
minDist = 4; %in pixels
scaleBar = 2; %in um
ext = '.omeTif';
path2File= 'E:\Data\Leuven Data\2018\ZHao\TestCode\400 nm AuNPs 1064 nm laser stepwise - 4_1';
%% Create the Movie object
switch ext
    case '.mp4'
        [file2Analyze] = Core.Movie.getFileInPath(path2File, ext);
        p2file = [file2Analyze.folder filesep file2Analyze.name];
        v = VideoReader(p2file);
        nFrames = floor(v.Duration*v.FrameRate);
        
        tic
        fullStackIn = zeros(v.Height,v.Width,nFrames);
        for i = 1:nFrames
            frame = readFrame(v);
            
            fullStackIn(:,:,i) = rgb2gray(frame);
        end
        toc
         
    otherwise
        info.type = 'transmission';
        myMov = Core.Movie(path2File,info,'.mp4');
        myMov.giveInfo;
        
        myMov.cropIm;

        fullStack = myMov.getFrame;

        fullStackIn = imcomplement(fullStack.Cam1);

end

%% detection of the center of the beads
frame =1;
%get the domaine
%[pos] = goldProj.simpleRegMaxDetection (fullStackIn(:,:,1),nParticles,minDist);
[pos] = Localization.smDetection(double(fullStackIn(:,:,frame)),minDist,2,12);
pos = round(pos);
im = double(fullStackIn(:,:,1));
[idx] = sub2ind(size(fullStackIn(:,:,1)),pos(:,1),pos(:,2));
bright = im(idx);

[~,idx2] = maxk(bright,nParticles);
[row,col] = ind2sub(size(fullStackIn(:,:,1)),idx(idx2));

x0 = col;
y0 = row;

pos = [y0 x0];
cropPos = round(mean(pos));

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

for i = 1:nFrames
    currentFrame = double(fullStackIn(:,:,i));
    
    %initial detection
%    [pos] = goldProj.simpleRegMaxDetection (currentFrame,nParticles,minDist);
    [pos] = Localization.smDetection(double(currentFrame),minDist,4,24);
    pos = round(pos);
    im = double(fullStackIn(:,:,1));
    [idx] = sub2ind(size(fullStackIn(:,:,1)),pos(:,1),pos(:,2));
    bright = im(idx);

    [~,idx2] = maxk(bright,nParticles);
    [row,col] = ind2sub(size(fullStackIn(:,:,1)),idx(idx2));

   
    x0 = col;
    y0 = row;
% %     x0 = pos(:,2);
%     y0 = pos(:,1);
    
    [gPar,resnorm,res] = Localization.Gauss.MultipleFitting(currentFrame,x0,y0,dom,nParticles); 
    F = Localization.Gauss.MultipleGauss(gPar, dom,nParticles);
    
    if i>1
        
        newOrder = goldProj.simpleTracking(gPar(5:end),prevPos);
        %reshaping and sorting with new order
        gPos = reshape(gPar(5:end),[2,nParticles]);
        gPos = gPos(:,newOrder);
        prevPos = reshape(gPos,[1,nParticles*2]);
        gPos = reshape(gPos,[1,2,nParticles]);
        data2Store(i,:,:) = gPos;
        clear gPos
        
    else
        
        gPos = reshape(gPar(5:end),[2,nParticles]);
        gPos = reshape(gPos,[1,2,nParticles]);
        data2Store(i,:,:) = gPos;
        %store the gPar in the prev for checking particle order
        prevPos = gPar(5:end);
        clear gPos
    end
    
    fitMov(:,:,i) = F;
    waitbar(i/nFrames,h,'Fitting Data')
end

filename = [path2File filesep 'LocalizationData.mat'];
save(filename,'data2Store');

close(h);
%% display figure
data2plot = data2Store *pxSize;
cm = [mean(mean(data2plot(:,1,:))) mean(mean(data2plot(:,2,:)))];
Fig = figure;
hold on
for i = 1 : nParticles
    scatter(data2plot(:,1,i)-cm(1),data2plot(:,2,i)-cm(2),20,'filled')
    
end
axis image
xlabel('X Position (px, 1px = 95nm)')
ylabel('Y Position (px, 1px = 95nm)')
title('All localized spot');
filename = [path2File filesep 'LocalizationDensity.fig'];
saveas(Fig,filename);
%% Cropping full stack
switch ext
    case '.mp4'
        fullStack = fullStackIn;
    otherwise
        fullStack = fullStack.Cam1(cropPos(1)-delta:cropPos(1)+delta, cropPos(2)-delta:cropPos(2)+delta,:);
end
%% MovieMaker

frameRate = 30;
filename = [path2File filesep 'TrackMovie.gif'];

goldProj.makeTraceMovie(data2Store,fullStack,filename,frameRate,scaleBar);

%% Calculate correlation

%xcorr??

RXY = corrcoef(x,y);

%process MP4;
