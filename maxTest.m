clc
clear 
close all
%% User input
path2File = 'D:\Documents\Unif\PhD\2018-Data\03 - March\MaxData\ECM remodeling\OnlyIntensity\testMax.tif';
threshold = 0.2;
%% Loading

info = Load.Movie.tif.getinfo(path2File);
frame = 1:info.Frame_n;

IM = Load.Movie.tif.getframes(path2File,frame);

figure(1)
title('Detecting cell countour')
subplot(2,4,1)
imagesc(IM(:,:,1))
title('Raw Data')
colormap('hot')
axis image
%% Segmenting cell

% size of gauss filter
S = 1;
% size of pixel in z vs x/y
pixZ  = 4;
zFactor = 2;
sigma = [S,S,S*zFactor/pixZ];
IMs = imgaussfilt3(IM, sigma);
disp('DONE with filtering ------------')


BW = imbinarize(IM,threshold);

figure(1)
subplot(2,4,2)
imagesc(BW(:,:,1))
title('Binary image')
colormap('hot')
axis image
contour = cell(size(BW,3),1);

for i = 1 : size(BW,3)
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
    %only keep largest countour
    if isempty(idx2OuterContour)
        cCountour = [];
    else
        cContour = cContour{idx2OuterContour};
    end
    contour{i} = cContour;
end

%plotting
figure(1)
subplot(2,4,3)
imagesc(BW(:,:,1))
hold on
colormap('hot')
axis image
title('Cleaned binary image')

figure(1)
subplot(2,4,4)
imagesc(IM(:,:,1))
hold on
plot(countour{1}(:,2),countour{1}(:,1),'w','LineWidth',3)
colormap('hot')
axis image
title('Cell Countour')

figure(2)
hold on
for i=1:size(IMs,3)
    if ~isempty(countour{i})
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


