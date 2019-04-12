clc
clear 
close all
%% User input
path2File = 'D:\Documents\Unif\PhD\2018-Data\03 - March\MaxData\ECM remodeling\OnlyIntensity\testMax.tif';
thresh1 = 0.20;
thresh2 = 0.20;
exF = 40;

%
%
%  A = surfaceArea(shp)
%
%

%% Loading

info = Load.Movie.tif.getinfo(path2File);
frame = 1:info.Frame_n;

IM = Load.Movie.tif.getframes(path2File,frame);
%% pre-processing
% size of gauss filter
S = 1;
% size of pixel in z vs x/y
pixZ  = 4;
zFactor = 2;
sigma = [S,S,S*zFactor/pixZ];
IMs = imgaussfilt3(IM, sigma);
disp('DONE with filtering ------------')

%% 3D segmentation

[gBW,aBW] = imSegmentation.segmentStack(IMs,'threshold',thresh1);
disp('DONE with segmenting ------------') 
%% test segmenting cell
  
   doplot = true; 
   [contour] = maxProj.getCellContour(gBW(:,:,exF),IM(:,:,exF),thresh1,doplot);
%% Segmenting cell
Contour = cell(size(IM,3),1);
outContour = cell(size(IM,3),1);
for i = 1 : size(IM,3)
        
       cIM = IM(:,:,i);
       cBW = gBW(:,:,i);
       doplot = false; 
      [contour,outCont] = maxProj.getCellContour(cBW,cIM,thresh1,doplot);
      Contour{i} = contour;
      outContour{i} = outCont;
     
end
%%
figure(2)
hold on
for i=1:size(IMs,3)
    if ~isempty(Contour{i})
        plot3(Contour{i}(:,2),Contour{i}(:,1),repmat(i,1,length(Contour{i}(:,1))),'k')
    end
end
%% inpsect contour
% figure(4)
% 
% for i=1:size(IMs,3)
%     hold on
%     imagesc(IMs(:,:,i))
%     if ~isempty(Contour{i})
%         plot(outContour{i}(:,2),outContour{i}(:,1),'r')
%         plot(Contour{i}(:,2),Contour{i}(:,1),'w')
%     end
%     axis image
%     i
%     waitforbuttonpress
% end
%% cleaning Contour
fContour = [];
for z = 1 :size(IMs,3)
    if ~isempty(Contour{z})
        fContour = [fContour ; Contour{z}(:,1) Contour{z}(:,2) ones(length(Contour{z}(:,1)),1)*z];
    end
end
%% get 3D contour
% test = scatteredInterpolant(fContour(:,2),fContour(:,2),fContour(:,1),ones(length(fContour(:,1)),1));
% [xq,yq,zq] = meshgrid(1:size(IMs,1),1:size(IMs,2),1:size(IMs,3));
% vq = test(xq,yq,zq);
% 
% 
% bd = boundary(fContour,1);
% interpCont(:,1) =fContour(bd(:,1),1);  
% interpCont(:,2) =fContour(bd(:,2),2);
% interpCont(:,3) = fContour(bd(:,3),3);

%% building 3D mask
mask3D = zeros(size(IMs));
for i = 1:size(IMs,3)
    idx = find(fContour(:,3)==i);
    mask3D(:,:,i) = poly2mask(fContour(idx,2),fContour(idx,1),size(IM,1),size(IM,2));
    
end

figure
subplot(1,2,1)
imagesc(squeeze(mean(mask3D,1)))
title('XZ Crosssection')

subplot(1,2,2)
imagesc(squeeze(mean(mask3D,3)))
title('XY Crosssection')
axis image
%% cleaninng mask
se = strel('cuboid',[1,1,2]);
%se = offsetstrel('ball',1,2);
mask3D =imdilate(mask3D,se);

mask3D = imerode(mask3D,se);
figure
subplot(1,2,1)
imagesc(squeeze(mean(mask3D,1)))
title('XZ Crosssection')

subplot(1,2,2)
imagesc(squeeze(mean(mask3D,3)))
title('XY Crosssection')
axis image
%% further cleaning
maxArea = regionprops(mask3D(:,:,40),'area');
thresh = 0.2;
for i = 1:size(IMs,3)
    cMask = mask3D(:,:,i);
    cArea = regionprops(cMask,'area');
    if ~isempty(cArea)
        if cArea.Area <thresh*maxArea.Area 
            mask3D(:,:,i) = zeros(size(mask3D,1),size(mask3D,2));
        end
    end
end

figure
subplot(1,2,1)
imagesc(squeeze(mean(mask3D,1)))

title('XZ Crosssection')

subplot(1,2,2)
imagesc(squeeze(mean(mask3D,3)))
title('XY Crosssection')
axis image
%% Segmenting the rest of the image 
% 
% mask = poly2mask(cContour(:,2),cContour(:,1),size(IM,1),size(IM,2));
% mask = ~mask;
% figure(1)
% subplot(2,4,5)
% imagesc(mask)
% colormap('hot');
% axis image;
% title('Mask');
threshold = 0.6;
diskDim = 4;
[gBW,aBW] = imSegmentation.segmentStack(IMs,'threshold',threshold,'diskDim',diskDim);
invMask3D = ~mask3D;
BWadapt  = aBW.*invMask3D;
BWglobal = gBW.*invMask3D;

% figure(1)
% subplot(2,4,6)
% imagesc(BWadapt(:,:,40))
% colormap('hot');
% axis image;
% title('Binary pores');

%% calculate the intensity-distance relationship

%pMask = poly2mask(cContour(:,2),cContour(:,1),size(IM,1),size(IM,2));
frameLimit1 = [ones(size(IMs,1),1),(1:size(IMs,2))'];
frameLimit2 = [ones(size(IMs,1),1)*size(IMs,1),(1:size(IMs,2))'];
frameLimit3 = [(1:size(IMs,1))',ones(size(IMs,2),1)];
frameLimit4 = [(1:size(IMs,1))',ones(size(IMs,2),1)*size(IMs,2)];
frameLimit  = [frameLimit1;frameLimit2;frameLimit3;frameLimit4];
idx2Lim = sub2ind(size(IMs),frameLimit(:,1),frameLimit(:,2));

%szDisk = 10;
cROI = mask3D;
pMask = mask3D;
tmp = mean(pMask,3);
tmp(tmp>0) = 1;
cont = bwboundaries(tmp);
cBound = cont{1};

%calculate the limit of the cell in X-Y
minX = min(cBound(:,2));
maxX = max(cBound(:,2));
minY = min(cBound(:,1));
maxY = max(cBound(:,1));

px = 10;
%SE = offsetstrel('ball',px,2*px);
SE = strel('cuboid',[px,px,px/2]);
pNPx = 1;
cNPx = 2;
count = 1;
step = [];
polDensity = [];
polIntensity = [];
fROI = ones(size(mask3D));
tROI = ones(size(mask3D));
medPolIntensity = [];
avgInt = 1;
figure
hold on
while ~all(cROI(:)==0)
    if length(step)<2
        pNPx =0;
    else
        pNPx = sum(sum(pROI));
    end

    cMask = imdilate(pMask,SE);
    cROI = cMask - pMask;
    if ~all(cROI(:)==0)
    cNPx = sum(sum(sum(cROI)));
    polDensity = [polDensity sum(sum(sum(double(BWadapt).*cROI)))/cNPx];
    avgInt = sum(sum(sum(double(IMs).*cROI)))/cNPx;
    polIntensity = [polIntensity avgInt];
    medPolIntensity = [medPolIntensity median(median(median(double(IMs).*cROI)))];
       
    tmp = mean(cROI,3);
    tmp(tmp>0) = 1;
    cont = bwboundaries(tmp);
    cBound = cont{1};

    %calculate the limit of the cell in X-Y
    cMinX = min(cBound(:,2));
    cMaxX = max(cBound(:,2));
    cMinY = min(cBound(:,1));
    cMaxY = max(cBound(:,1));

    DX = round(mean([minX-cMinX cMaxX-maxX]));
    DY = round(mean([minY-cMinY cMaxY-maxY]));

    step(count) = round(mean([DX(DX>0) DY(DY>0)]));
  
    fROI(cROI==1) = avgInt;
    pMask = cMask;    
    pROI = cROI;
    count = count+1;
    
    tROI(cROI==1) = count;
    imagesc(tROI(:,:,40));
   
    end
    
end
%%
figure(3)

%  subplot(1,2,1)
%  plot(step,polDensity);
%  title('Mak');
%  xlabel('Distance from cell - Pixel')
%  ylabel('Proportion of black/white pixel')
% y = ones(size(step))*median(median(median(double(IMs))));
% y2 = ones(size(step))*mean(mean(mean(double(IMs))));
subplot(1,2,2)
plot(step,polIntensity);
xlabel('Distance from cell - Pixel');
ylabel('Average intensity per pixel');
hold on
plot(step,y)
plot(step,y2)
legend({'Avg Intensity per pixel','Median intensity full image', 'Mean intensity full image',})


% %% calculate the intensity-distance relationship
% 
% pMask = poly2mask(cContour(:,2),cContour(:,1),size(IM,1),size(IM,2));
% frameLimit1 = [ones(size(IMs,1),1),(1:size(IMs,2))'];
% frameLimit2 = [ones(size(IMs,1),1)*size(IMs,1),(1:size(IMs,2))'];
% frameLimit3 = [(1:size(IMs,1))',ones(size(IMs,2),1)];
% frameLimit4 = [(1:size(IMs,1))',ones(size(IMs,2),1)*size(IMs,2)];
% frameLimit  = [frameLimit1;frameLimit2;frameLimit3;frameLimit4];
% idx2Lim = sub2ind(size(IMs),frameLimit(:,1),frameLimit(:,2));
% 
% szDisk = 10;
% cROI = pMask;
% SE = strel('disk',szDisk);
% pNPx = 1;
% cNPx = 2;
% count = 1;
% step = [];
% polDensity = [];
% polIntensity = [];
% 
% while all(ismember(find(cROI),idx2Lim)==0)
%     if length(step)<2
%         pNPx =0;
%     else
%         pNPx = sum(sum(pROI));
%     end
%         
%     cMask = imdilate(pMask,SE);
%     cROI = cMask - pMask;
%     figure(2)
%     hold on
%     imagesc(cROI)
%     colormap('hot');
%     axis image
%     drawnow
%     pause(0.1);
%     hold off
%     
%     cNPx = sum(sum(cROI));
%     polDensity = [polDensity sum(sum(double(BWadapt).*cROI))/cNPx];
%     polIntensity = [polIntensity sum(sum(double(IMs).*cROI))/cNPx];
%     
%     step(count) = count*(szDisk-1);
%     
%     pMask = cMask;    
%     pROI = cROI;
%     count = count+1;
% end
% 
% figure(3)
% 
% subplot(1,2,1)
% plot(step,polDensity);
% 
% y = ones(size(step))*median(median(double(IMs)));
% y2 = ones(size(step))*mean(mean(double(IMs)));
% subplot(1,2,2)
% plot(step,polIntensity);
% hold on
% plot(step,y)
% plot(step,y2)
% 

