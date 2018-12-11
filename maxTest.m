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

figure

imagesc(IM)
colormap('hot')
%% Segmenting cell

BW = imbinarize(IM,threshold);

figure

imagesc(BW)
colormap('hot')

BWL = bwlabel(BW);

test = regionprops(BW,'Area');
[biggestArea,idx] = max(cell2mat({test.Area}));

BW(BWL~=idx) = 0;

% Clean up boundary
se = strel('disk',10);
BW = imclose(BW,se);


figure(3)
imagesc(BW)
hold on
colormap('hot')

contour = bwboundaries(BW);
figure(4)
imagesc(IM)
hold on
plot(contour{1}(:,2),contour{1}(:,1),'w','LineWidth',3)
colormap('hot')

%% Segmenting the rest of the image 

%% calculate the ellipticity-distance relationship


