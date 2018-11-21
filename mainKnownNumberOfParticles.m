clear 
clc
close all
%% User input
delta = 40;
nParticles = 4;

%% Create the Movie object

path2File= 'E:\Data\Leuven Data\2018\ZHao\TestCode\400 nm AuNPs 1064 nm laser stepwise - 4_1';
info.type = 'transmission';
myMov = Core.Movie(path2File,info);
myMov.giveInfo;

%%
myMov.cropIm;

%% get the data from the movie
fullStack = myMov.getFrame;

%% inversion of the scale cam1 extraction

fullStackIn = imcomplement(fullStack.Cam1);
%% detection of the center of the beads
%get the domaine
[pos] = goldProj.simpleRegMaxDetection (fullStackIn(:,:,1),nParticles);
cropPos = round(mean(pos));

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
    [pos] = goldProj.simpleRegMaxDetection (currentFrame,nParticles);
    x0 = pos(:,2);
    y0 = pos(:,1);
    
    [gPar,resnorm,res] = Localization.Gauss.MultipleFitting(currentFrame,x0,y0,dom,nParticles); 
    F = Localization.Gauss.MultipleGauss(gPar, dom,nParticles);
    
    %let us check if particle order somehow flip
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

Fig = figure;
hold on
for i = 1 : nParticles
    scatter(data2Store(:,1,i),data2Store(:,2,i),20,'filled')
    
end
axis image
xlabel('X Position (px, 1px = 95nm)')
ylabel('Y Position (px, 1px = 95nm)')
title('All localized spot');
filename = [path2File filesep 'LocalizationDensity.fig'];
saveas(Fig,filename);
%% Cropping full stack
fullStack = fullStack.Cam1(cropPos(1)-delta:cropPos(1)+delta, cropPos(2)-delta:cropPos(2)+delta,:);
%% MovieMaker

frameRate = 30;
filename = [path2File filesep 'TrackMovie.gif'];

goldProj.makeTraceMovie(data2Store,fullStack,filename,frameRate);


