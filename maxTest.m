clc
clear 
close all
%% User input
path2File = 'D:\Documents\Unif\PhD\2018-Data\03 - March\MaxData\ECM remodeling\OnlyIntensity\testMax.tif';
threshold = 0.20;
exF = 20;
%% Loading

info = Load.Movie.tif.getinfo(path2File);
frame = 1:info.Frame_n;

IM = Load.Movie.tif.getframes(path2File,frame);

figure(1)
title('Detecting cell countour')
subplot(2,4,1)
imagesc(IM(:,:,exF))
title('Raw Data')
colormap('hot')
axis image
%% test segmenting cell
    threshold = 0.15;
    BW(:,:,20) = imbinarize(IM(:,:,20),threshold);
    cBW = BW(:,:,20);
    cBWL = bwlabel(cBW);

    %Get the largest area
    cBWarea = regionprops(cBW,'Area');
    [~,idx2BiggestArea] = max(cell2mat({cBWarea.Area}));
    if isempty(idx2BiggestArea)
    else
    %kill all the other area found
        cBW(cBWL~=idx2BiggestArea) = 0;
    end
    figure 
    imagesc(cBW)
    colormap('hot');
    % Clean up boundary
    se = strel('disk',10);
    cBW = imclose(cBW,se);
    
    cContour = bwboundaries(cBW);
    
    %find largest contour
    contourSize = cellfun(@size,cContour,'UniformOutput',false);
    [tmp] = cellfun(@max,contourSize,'UniformOutput',false);
    [~,idx2OuterContour] = max(cell2mat(tmp));
    %only keep largest countour
    if isempty(idx2OuterContour)
        cContour = [];
    else
        cContour = cContour{idx2OuterContour};
    end
    
    mask = poly2mask(cContour(:,2),cContour(:,1),size(IM,1),size(IM,2));
    
    nIM = double(IM(:,:,20)).*mask;
    %%
    figure 
    imagesc(nIM)
    nIM(nIM==0) = min(nIM(nIM~=0));
    
    figure
    imagesc(nIM)
    
  %  [gBW,aBW] = imSegmentation.segmentStack(nIM,'Threshold',0.0001);
    Thresh = adaptthresh(nIM,0.1);
%%
threshold = 0.2;
%gBW = imbinarize(nIM,0.9);
gBW = nIM;

gBW(gBW<threshold*max(max(gBW))) = 0;
gBW(gBW>=threshold*max(max(gBW))) =1;
figure 
imagesc(gBW);
    
 %% cleaning
nBW = gBW;
se = strel('disk',3);
nBW = imclose(nBW,se);
figure 
imagesc(nBW)
nBW = ~nBW;
nBWL = bwlabel(nBW);
nBWarea = regionprops(nBW,'Area');
[val,idx2BiggestArea] = maxk(cell2mat({nBWarea.Area}),2);
if isempty(idx2BiggestArea)
else
%kill all the other area found
    [~,idx] = min(val);
    nBW(nBWL~=idx2BiggestArea(idx))= 0;
end
 
figure
imagesc(nBW)
% Clean up boundary
se = strel('disk',10);
nBW = imclose(nBW,se);
figure
imagesc(nBW)
cContour = bwboundaries(nBW);
%find largest contour
contourSize = cellfun(@size,cContour,'UniformOutput',false);
[tmp] = cellfun(@max,contourSize,'UniformOutput',false);
[~,idx2OuterContour] = max(cell2mat(tmp));
%only keep largest countour
if isempty(idx2OuterContour)
    cContour = [];
else
    cContour = cContour{idx2OuterContour};
end

figure
imagesc(IM(:,:,20))
hold on
plot(cContour(:,2),cContour(:,1),'w','LineWidth',3)

%% Segmenting cell

% % size of gauss filter
% S = 1;
% % size of pixel in z vs x/y
% pixZ  = 4;
% zFactor = 2;
% sigma = [S,S,S*zFactor/pixZ];
% IMs = imgaussfilt3(IM, sigma);
% disp('DONE with filtering ------------')
IMs = IM;

BW1 = imbinarize(IM,threshold);

figure(1)
subplot(2,4,2)
imagesc(BW1(:,:,exF))
title('Binary image')
colormap('hot')
axis image
contour = cell(size(IM,3),1);
BW = BW1;%copy for testing
for i = 1 : size(IM,3)
    BW(:,:,i) = imbinarize(IM(:,:,i),threshold);
    cBW = BW(:,:,i);
    cBWL = bwlabel(cBW);

    %Get the largest area
    cBWarea = regionprops(cBW,'Area');
    [~,idx2BiggestArea] = max(cell2mat({cBWarea.Area}));
    if isempty(idx2BiggestArea)
    else
    %kill all the other area found
        cBW(cBWL~=idx2BiggestArea) = 0;
    end
    % Clean up boundary
    se = strel('disk',10);
    cBW = imclose(cBW,se);
    
    cContour = bwboundaries(cBW);
    
    %find largest contour
    contourSize = cellfun(@size,cContour,'UniformOutput',false);
    [tmp] = cellfun(@max,contourSize,'UniformOutput',false);
    [~,idx2OuterContour] = max(cell2mat(tmp));
    [~,idx2InnerContour] = min(cell2mat(tmp));
    %only keep largest countour
    if isempty(idx2OuterContour)
        cContour = [];
    else
        cContour = cContour{idx2OuterContour};
    end
    contour{i} = cContour;
end

%plotting
figure(1)
subplot(2,4,3)
imagesc(BW(:,:,exF))
hold on
colormap('hot')
axis image
title('Cleaned binary image')

figure(1)
subplot(2,4,4)
imagesc(IM(:,:,exF))
hold on
plot(contour{exF}(:,2),contour{exF}(:,1),'w','LineWidth',3)
colormap('hot')
axis image
title('Cell Countour')

figure(2)
hold on
for i=1:size(IMs,3)
    if ~isempty(contour{i})
        plot3(contour{i}(:,2),contour{i}(:,1),repmat(i,1,length(contour{i}(:,1))),'k')
    end
end

%% Segmenting the rest of the image 

mask = poly2mask(cContour(:,2),cContour(:,1),size(IM,1),size(IM,2));
mask = ~mask;
figure(1)
subplot(2,4,5)
imagesc(mask)
colormap('hot');
axis image;
title('Mask');

disp('Now doing 3D gauss filtering this can take about 3 minutes')


diskDim = 4;
[gBW,aBW] = imSegmentation.segmentStack(IMs,'threshold',threshold,'diskDim',diskDim);

BWadapt  = aBW.*mask;
BWglobal = gBW.*mask;

figure(1)
subplot(2,4,6)
imagesc(BWadapt)
colormap('hot');
axis image;
title('Binary pores');

%% calculate the ellipticity-distance relationship

pMask = poly2mask(cContour(:,2),cContour(:,1),size(IM,1),size(IM,2));
frameLimit1 = [ones(size(IMs,1),1),(1:size(IMs,2))'];
frameLimit2 = [ones(size(IMs,1),1)*size(IMs,1),(1:size(IMs,2))'];
frameLimit3 = [(1:size(IMs,1))',ones(size(IMs,2),1)];
frameLimit4 = [(1:size(IMs,1))',ones(size(IMs,2),1)*size(IMs,2)];
frameLimit  = [frameLimit1;frameLimit2;frameLimit3;frameLimit4];
idx2Lim = sub2ind(size(IMs),frameLimit(:,1),frameLimit(:,2));

szDisk = 10;
cROI = pMask;
SE = strel('disk',szDisk);
pNPx = 1;
cNPx = 2;
count = 1;
step = [];
polDensity = [];
polIntensity = [];

while all(ismember(find(cROI),idx2Lim)==0)
    if length(step)<2
        pNPx =0;
    else
        pNPx = sum(sum(pROI));
    end
        
    cMask = imdilate(pMask,SE);
    cROI = cMask - pMask;
    figure(2)
    hold on
    imagesc(cROI)
    colormap('hot');
    axis image
    drawnow
    pause(0.1);
    hold off
    
    cNPx = sum(sum(cROI));
    polDensity = [polDensity sum(sum(double(BWadapt).*cROI))/cNPx];
    polIntensity = [polIntensity sum(sum(double(IMs).*cROI))/cNPx];
    
    step(count) = count*(szDisk-1);
    
    pMask = cMask;    
    pROI = cROI;
    count = count+1;
end

figure(3)

subplot(1,2,1)
plot(step,polDensity);

y = ones(size(step))*median(median(double(IMs)));
y2 = ones(size(step))*mean(mean(double(IMs)));
subplot(1,2,2)
plot(step,polIntensity);
hold on
plot(step,y)
plot(step,y2)


