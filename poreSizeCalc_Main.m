clear;
clc;
close all;
%% User Input
pxSize = 100; %in nm
Threshold = 0.5; %number between 0 and 1 (% of max intensity)
bigPores = 25; %in px ("draw a line" between small pores and big pores);
nBins = 20; %For Log histogram
nFrame =10; %n of frame to analyze
%% Loading Data
pxArea = pxSize*pxSize*1e-6; %in µm^2
bigPores = bigPores*pxArea;
mainFolderName = uigetdir;
assert(ischar(mainFolderName),'User canceled the selection of file, excecution aborted');
%extract the name of the current folder
idx = strfind(mainFolderName,'\') ;
currentFolderName = mainFolderName(idx(end)+1:end) ;
%Remove dots from the name
currentFolderName = regexprep(currentFolderName,'\.','_');

%Extract the part of the folder that is a tif file
Folder_Content = dir(mainFolderName);
index2Images   = contains({Folder_Content.name},'.tif');
images2Analyze = Folder_Content(index2Images);
%% Displaying First image
stack2Load = 1;
frame2Load = [10 50 100 140];

path2Stacks = strcat(images2Analyze(stack2Load).folder,'\');
p2file      = strcat(path2Stacks,images2Analyze(stack2Load).name);
warning('off');
fileInfo    = loadMovie.tif.getinfo(p2file);

for n=1:length(frame2Load)
IM     = loadMovie.tif.getframes(p2file, frame2Load(n));
warning('on');

IM = imgaussfilt(IM,2);
T = adaptthresh(IM,Threshold,'ForegroundPolarity','dark','Statistic','mean');

BW = imbinarize(IM,T);
BW = bwareaopen(BW,4);
se = strel('disk',2);
BW = imclose(BW,se);

IM2 = IM;
IM3 = IM;
IM2(BW) = 0;
IM3(~BW) = mean(mean(IM(~BW))); 
IM4 = IM2;
IM4(IM4>0) = 1;

H0 = figure(n);
hold(gca,'on')
subplot(1,3,1)

imagesc(IM)
axis image

subplot(1,3,2)
imagesc(IM2)
axis image

subplot(1,3,3)
imagesc(IM3)
axis image
title(sprintf('Frame %d',frame2Load(n)))
hold off

fileName0 = sprintf('%s%s%s-PoresFrame-%d_Threshold-%d',mainFolderName,'\',currentFolderName,frame2Load(n),Threshold*100);
savefig(H0,fileName0)

end

%% Looping through the Data

h = waitbar(0);
bins = zeros(nBins,size(images2Analyze,1)); %Store Bins
occurrences = bins; % Store occurences

for j = 1:size(images2Analyze,1)
    hMessage = sprintf('Loading image stack number %d/%d',j,size(images2Analyze,1));
    waitbar(0,h,hMessage);
    path2Stacks = strcat(images2Analyze(j).folder,'\');
    p2file      = strcat(path2Stacks,images2Analyze(j).name);
    warning('off');
    fileInfo    = loadMovie.tif.getinfo(p2file);
    %frames2Load = 1:1:fileInfo.Frame_n;
    frames2Load  = 1:1:nFrame; 
    IMStack     = loadMovie.tif.getframes(p2file, frames2Load);
    warning('on');
    Results = struct('numPores',cell(1,size(images2Analyze,1)),'numBigPores',cell(1,size(images2Analyze,1)),'Area',...
        cell(1,size(images2Analyze,1)),'areaBigPores',cell(1,size(images2Analyze,1)),...
        'meanArea',cell(1,size(images2Analyze,1)),'meanAreaBigPores',cell(1,size(images2Analyze,1)),...
        'medArea',cell(1,size(images2Analyze,1)),'medAreaBigPores',cell(1,size(images2Analyze,1)),...
        'CVArea',cell(1,size(images2Analyze,1)),'CVAreaBigPores',cell(1,size(images2Analyze,1)),...
        'MajorAxis',cell(1,size(images2Analyze,1)),'MinorAxis',cell(1,size(images2Analyze,1)),...
        'meanElip',cell(1,size(images2Analyze,1)),'medElip',cell(1,size(images2Analyze,1)));
    dataForHistogram = [];
    hMessage = sprintf('Analysis of Stack Number %d/%d',j,size(images2Analyze,1));
    for i=1:size(IMStack,3)
    waitbar(i/size(IMStack,3),h,hMessage);
        % Loading image number i
        IM = IMStack(:,:,i);
        I  = double(IM);
        % Normalize the image
        I  = I./max(max(I));
        % Gaussian filtering (smooth image)
        I  = imgaussfilt(I,2);
        % Get an adaptive threshold (~depends on local intensity levels)
        T  = adaptthresh(I,Threshold,'ForegroundPolarity','dark','Statistic','mean');
        % Binarization of image and cleaning
        BW = imbinarize(I,T);
        BW = bwareaopen(BW,4);
        se = strel('disk',2);
        BW = imclose(BW,se);
        BW_pores = ~BW;
        %Remove pores that are in touch with the edges(uncertainty about
        %size)
        BW_pores = imclearborder(BW_pores);
        
        %Get properties of the pores on the image.
        stats = regionprops(BW_pores,'Area','MajorAxisLength',...
        'MinorAxisLength');
        
        if isempty(stats) %avoid concatenation errors due to empty cell later
            clear stats;
            stats.Area = NaN;
            stats.MajorAxisLength = NaN;
            stats.MinorAxisLength = NaN;
            Results(i).numPores = 0;
            
        else
        Results(i).numPores = size(stats,1);
        Results(i).Area = [stats.Area].*pxArea;
        Results(i).MajorAxis = [stats.MajorAxisLength];
        Results(i).MinorAxis = [stats.MinorAxisLength];
        end
  
        Results(i).numBigPores = size(Results(i).Area(Results(i).Area>bigPores),2);
        Results(i).areaBigPores = Results(i).Area(Results(i).Area>bigPores);
        
        if isempty(Results(i).areaBigPores)%avoid concatenation errors due to empty cell later
            Results(i).areaBigPores = 0;
        end
        
        Results(i).meanArea = mean(Results(i).Area);
        Results(i).meanAreaBigPores = mean(Results(i).Area(Results(i).Area>bigPores));
        Results(i).medArea  = median(Results(i).Area);
        Results(i).medAreaBigPores = median(Results(i).Area(Results(i).Area>bigPores));
        Results(i).CVArea = std(Results(i).Area)/Results(i).meanArea;
        Results(i).CVAreaBigPores = std(Results(i).areaBigPores)/Results(i).meanAreaBigPores;
        Results(i).Ellipticity = [stats.MinorAxisLength]./[stats.MajorAxisLength];
        Results(i).meanElip = mean(Results(i).Ellipticity);
        Results(i).medElip = median(Results(i).Ellipticity);
        dataForHistogram = [dataForHistogram [stats.Area].*pxArea];
    end
   
    %% Save results to excel sheets

    fileNameExcel = sprintf('%s%sResultsSummary',mainFolderName,'\');
    matName =  regexprep(images2Analyze(j).name,'\.','_');
    fileNameMat   = sprintf('%s%s%s-FullResults',mainFolderName,'\',...
        matName);
    save(fileNameMat,'Results');
    sheetName = sprintf('Stack%d',j);
    Title=[ 'numPores    ';'numBigPores ';'meanArea    ';'meanBigPores';...
            'medArea     ';'medBigPores ';'CVArea      ';'CVAreaBig   ';...
            'meanEllip   ';'medEllip    '];
        Titletable=cellstr(Title)';
        xlswrite(fileNameExcel,Titletable,sheetName,'A1:H1');

    Result2print = [[Results.numPores]', [Results.numBigPores]',...
        [Results.meanArea]',[Results.meanAreaBigPores]',[Results.medArea]',...
        [Results.medAreaBigPores]',[Results.CVArea]',[Results.CVAreaBigPores]',...
        [Results.meanElip]',[Results.medElip]'];
    
        sizeRes = size(Results,2);
        Range = sprintf('A2:H%d',sizeRes+1);
        
        xlswrite(fileNameExcel,Result2print, sheetName,Range);
    
    [midBin, ~,occurrence]=Misc.lnbin(dataForHistogram,20);
    bins(:,j) = midBin;
    occurrences(:,j) = occurrence;
    histData(j).bins = midBin;
    histData(j).occurrences = occurrences;
end
 close(h);
h = msgbox('The Data were succesfully saved !', 'Success');
%% Plotting
H1 = figure(10);
hold (gca,'on')
title(currentFolderName)
set(gca,'YScale','log');
set(gca,'XScale','log');
ylabel('Number of Pores');
xlabel('Pore area');
leg = cell(1,size(images2Analyze,1));
for i = 1:size(images2Analyze,1)
    plot(bins(:,i),occurrences(:,i));
    leg{i} = sprintf('Stack %d',i);

end

legend(leg)
hold (gca,'off')


H2 = figure(11);
hold(gca,'on')
title(sprintf('%s - Overview',currentFolderName))
scatter(median(bins,2),median(occurrences,2))
errorbar(median(bins,2),median(occurrences,2),std(occurrences,1,2),'LineStyle','none')
set(gca,'YScale','log');
set(gca,'XScale','log');
ylabel('Number of Pores');
xlabel('Pore area');
legend('median','Standard deviation')

hold(gca,'off')

histData(1).medBins = median(bins,2);
histData(1).medOcc  = median(occurrences,2);
histData(1).STD     = std(occurrences,1,2);
%% Saving figures & Data


fileName1 = sprintf('%s%s%s-AllCurves',mainFolderName,'\',currentFolderName);
savefig(H1,fileName1)

fileName2 = sprintf('%s%s%s-AverageCurve',mainFolderName,'\',currentFolderName);
savefig(H2,fileName2)

fileNameMat = sprintf('%s%s%s-histData',mainFolderName,'\',currentFolderName);
save(fileNameMat,'histData');

