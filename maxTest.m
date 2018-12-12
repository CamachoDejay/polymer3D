clc
clear 
close all
%% User input
path2File = 'C:\Users\Boris\Documents\MATLAB\data\max\3nm_3D_Series001_z28_ch00.tif';
threshold = 0.5;
%% Loading

info = Load.Movie.tif.getinfo(path2File);
frame = 1:info.Frame_n;

IM = Load.Movie.tif.getframes(path2File,frame);

figure(1)
title('Detecting cell countour')
subplot(2,4,1)
imagesc(IM)
title('Raw Data')
colormap('hot')
axis image
%% Segmenting cell

BW = imbinarize(IM,threshold);

figure(1)
subplot(2,4,2)
imagesc(BW)
title('Binary image')
colormap('hot')
axis image

BWL = bwlabel(BW);

%Get the largest area
BWarea = regionprops(BW,'Area');
[~,idx2BiggestArea] = max(cell2mat({BWarea.Area}));

%kill all the other area found
BW(BWL~=idx2BiggestArea) = 0;

% Clean up boundary
se = strel('disk',10);
BW = imclose(BW,se);

%plotting
figure(1)
subplot(2,4,3)
imagesc(BW)
hold on
colormap('hot')
axis image
title('Cleaned binary image')
contour = bwboundaries(BW);

%find largest contour
contourSize = cellfun(@size,contour,'UniformOutput',false);
[tmp] = cellfun(@max,contourSize,'UniformOutput',false);
[~,idx2OuterContour] = max(cell2mat(tmp));
%only keep largest countour
contour = contour{idx2OuterContour};

figure(1)
subplot(2,4,4)
imagesc(IM)
hold on
plot(contour(:,2),contour(:,1),'w','LineWidth',3)
colormap('hot')
axis image
title('Cell Countour')
%% Segmenting the rest of the image 

mask = poly2mask(contour(:,2),contour(:,1),size(IM,1),size(IM,2));
mask = ~mask;
figure(1)
subplot(2,4,5)
imagesc(mask)
colormap('hot');
axis image;
title('Mask');

disp('Now doing 3D gauss filtering this can take about 3 minutes')

% size of gauss filter
S = 1;
% size of pixel in z vs x/y
pixZ  = 4;
zFactor = 2;
sigma = [S,S,S*zFactor/pixZ];
IMs = imgaussfilt3(IM, sigma);
disp('DONE with filtering ------------')

diskDim = 4;
[gBW,aBW] = imSegmentation.segmentStack(IMs,'threshold',threshold,'diskDim',diskDim);

BWadapt  = aBW.*mask;
BWglobal = gBW.*mask;

figure(1)
subplot(2,4,6)
imagesc(BWadapt)
colormap('hot');
axis image;
title('Mask');

%% calculate the ellipticity-distance relationship


