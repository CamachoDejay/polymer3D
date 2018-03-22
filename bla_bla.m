%% Loading image 1
p2file = '/Users/rafa/Documents/MATLAB/data/test-pic.tif';
fileInfo = loadMovie.tif.getinfo(p2file);
 
IM = loadMovie.tif.getframes(p2file, 10);
 
%%
T = adaptthresh(IM,0.4,'ForegroundPolarity','dark','Statistic','mean');
 
BW = imbinarize(IM,T);
BW = bwareaopen(BW,4);
se = strel('disk',2);
BW = imclose(BW,se);

IM2 = IM;
IM3 = IM;
IM2(BW) = 0;
IM3(~BW) = mean(mean(IM(~BW)));
% IM3(BW)  = mean(mean(IM(BW)));

%%
figure(1)
subplot(1,2,1)
imagesc(IM)
axis image
a = gca;
a.XTickLabel = [];
a.YTickLabel = [];
 
subplot(1,2,2)
imagesc(BW)
axis image
a = gca;
a = gca;
a.XTickLabel = [];
a.YTickLabel = [];
% 
% subplot(1,3,3)
% imagesc(IM3)
% axis image

%%
figure(2)
% subplot(1,2,1)
imagesc(IM(1:200,1:200))
axis image
a = gca;
a.XTickLabel = [];
a.YTickLabel = [];


%% Loading image 1
p2file = '/Users/rafa/Documents/MATLAB/data/test-pic.tif';
fileInfo = loadMovie.tif.getinfo(p2file);
IM = loadMovie.tif.getframes(p2file, 10);

figure(1)
imagesc(IM);
axis image

[cropIM] = Misc.getImArea(IM);

figure(1)
clf
imagesc(cropIM);
axis image
%%
IM2proc = imgaussfilt(IM,3);
BWtest = imbinarize(IM2proc);
holes = ~BWtest;
holes = bwareaopen(holes,9);
figure(1)
subplot(1,2,1)
imagesc(IM)
axis image

subplot(1,2,2)
imagesc(holes)
axis image


%%

[L,n] = bwlabel(holes);

outData = zeros(n,3);
for i = 1:n
    tmpBW = L==i;
    [ fWidth ] = SDcalc.fastWidthBW( tmpBW );
    tmpSize = sum(sum(tmpBW));
    outData(i,1) = fWidth;
    outData(i,2) = tmpSize;
    
    [B] = bwboundaries(tmpBW,'noholes');
    Blength = cellfun(@length,B);
    [~, idx] = maxk(Blength,2);
    B = B(idx);
    assert(length(B)==1,'problems');
    boundary = B{1};
    [ vals, names ] = SDcalc.solidity( boundary' );
    outData(i,3) = vals{1};
    
end

%%

[B,L] = bwboundaries(holes,'noholes');

Blength = cellfun(@length,B);
% B(Blength<3) = [];

figure(1)
clf
imagesc(holes)
axis image

hold on
outSD = zeros(length(B),6);
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
   [~,~,a] = SDcalc.centroid_by_area( boundary' );
   outSD(k,1) = abs(a);
   [ vals, names ] = SDcalc.solidity( boundary' );
   outSD(k,2) = vals{1};
   [ fWidth ] = SDcalc.fastWidth( boundary' );
   outSD(k,3) = fWidth;
%    [ vals, names ] = SDcalc.contour_aspect_ratio_bb( boundary' );
%    outSD(k,4) = vals{1};
%    outSD(k,5) = vals{2};
%    outSD(k,6) = vals{3};
end
hold off
axis image