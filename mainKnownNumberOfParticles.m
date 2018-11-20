clear 
clc
close all
%% User input
chi2 = 80;
FWHM_pix = 4;
delta = 20;
nParticles = 2;

%% Create the Movie object

path2File= 'E:\Data\Leuven Data\2018\ZHao\TestCode\400 nm AuNPs 1064 nm laser stepwise - 2 polarization_1';
info.type = 'transmission';
myMov = Core.Movie(path2File,info);
myMov.giveInfo;

%%
myMov.showFrame(1,2);

%% get the data from the movie
fullStack = myMov.getFrame;

%% inversion of the scale cam1 extraction

fullStackIn = imcomplement(fullStack.Cam1);
%% detection of the center of the beads
test = double(fullStackIn(:,:,1)); %Conversion to double for detection
%get the domaine
x = 1:size(fullStackIn,2);
y = 1:size(fullStackIn,1);
[domX,domY] = meshgrid(x,y);

%rough estimate of bkg:
bkg = median(median(test));
test = test-bkg; %bkg subtraction
test(test<0.1*max(max(test))) = 0;
out = imregionalmax(test);

x0 = domX(out);
x0 = x0(x0~=0);
y0 = domY(out);
y0 = y0(y0~=0);
bright = test(out);
bright = bright(bright~=0);
%find the brightess maxima according to the number of particles specified
[~,idx] = maxk(bright,nParticles);
x0 = x0(idx);
y0 = y0(idx);
%[ pos, ~, ~] = Localization.smDetection(firstFrame, delta, FWHM_pix, chi2 );
pos = round(mean([y0,x0]));

%% Cropping Movie

fullStackIn = fullStackIn(pos(1)-delta:pos(1)+delta, pos(2)-delta:pos(2)+delta,:);
%% Data manipulation
nFrames = size(fullStackIn,3);
x = 1:size(fullStackIn,2);
y = 1:size(fullStackIn,1);
[domX,domY] = meshgrid(x,y);
dom(:,:,1) = domX;
dom(:,:,2) = domY;

data2Store = table(zeros(nFrames,1),zeros(nFrames,1),zeros(nFrames,1),zeros(nFrames,1),zeros(nFrames,1),...
    'VariableNames',{'Amp','x1','y1','x2','y2'});
fitMov = zeros(size(fullStackIn));
h = waitbar(0,'Fitting Data');
for i = 1:nFrames
    currentFrame = double(fullStackIn(:,:,i));
    
    %initial detection
    test = currentFrame;
    test  = test-bkg; %bkg subtraction
    test(test<0.1*max(max(test))) = 0;
    out = imregionalmax(test);
    x0 = domX(out);
    x0 = x0(x0~=0);
    y0 = domY(out);
    y0 = y0(y0~=0);
    bright = test(out);
    bright = bright(bright~=0);
   
    %find the brightess maxima according to the number of particles specified
    [~,idx] = maxk(bright,nParticles);
    x0 = x0(idx);
    y0 = y0(idx);
    
    [gPar,resnorm,res] = Localization.Gauss.MultipleFitting(currentFrame,x0,y0,dom,nParticles); 
    F = Localization.Gauss.MultipleGauss(gPar, dom,nParticles);
    
    %let us check if particle order somehow flip
    if i>1
        
        cX = [gPar(5) gPar(7)];
        cY = [gPar(6) gPar(8)];
        euclDist = sqrt((cX-prevX1).^2 + (cY-prevY1).^2);
        [~,idx] = min(euclDist);
        
        if idx == 1
            
            data2Store.x1(i)  = gPar(5);
            data2Store.y1(i)  = gPar(6);
            data2Store.x2(i)  = gPar(7);
            data2Store.y2(i)  = gPar(8);
            prevX1 = gPar(5);
            prevY1 = gPar(6);
            
        elseif idx == 2
            
            data2Store.x1(i)  = gPar(7);
            data2Store.y1(i)  = gPar(8);
            data2Store.x2(i)  = gPar(5);
            data2Store.y2(i)  = gPar(6);
            prevX1 = gPar(7);
            prevY1 = gPar(8);
        else
            warning('Could not decide which particle is which, stopping trace');
        end
    else
        data2Store.x1(i)  = gPar(5);
        data2Store.y1(i)  = gPar(6);
        data2Store.x2(i)  = gPar(7);
        data2Store.y2(i)  = gPar(8);
        prevX1 = gPar(5);
        prevY1 = gPar(6);
 
    end
    
    data2Store.Amp(i) = gPar(1);
    fitMov(:,:,i) = F;
    waitbar(i/nFrames,h,'Fitting Data')
end
close(h);

figure
scatter(data2Store.x1,data2Store.y1,20,'filled')
hold on
scatter(data2Store.x2,data2Store.y2,20,'filled')

%% Cropping full stack
fullStack = fullStack.Cam1(pos(1)-delta:pos(1)+delta, pos(2)-delta:pos(2)+delta,:);
%% MovieMaker

frameRate = 30;
filename = [path2File filesep 'TrackMovie.gif'];
Fig = figure;

for i = 1:nFrames
    
    cFrame = fullStack(:,:,i);

    imagesc(cFrame)
    hold on
    colormap('gray')
    axis image;
    
    plot(data2Store.x1(i),data2Store.y1(i),'b+')
    plot(data2Store.x1(1:i),data2Store.y1(1:i),'b-')
    
    plot(data2Store.x2(i),data2Store.y2(i),'r+')
    plot(data2Store.x2(1:i),data2Store.y2(1:i),'r-')
    
    
    set(gca,'visible','off');
    set(gcf,'color','w');
    drawnow;
    
    frame = getframe(Fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if i == 1

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);

    else

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');

    end
    
    hold off
    clf;
    
    
end


