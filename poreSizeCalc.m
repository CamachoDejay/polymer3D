clear;
clc;
close all;

%% User Input
pxSize = 100; %in nm
Threshold = 0.5; %number between 0 and 1 (% of max intensity)
bigPores = 25; %in px ("draw a line" between small pores and big pores);
nBins = 20; %For Log histogram
nFrame = 10; %n of frame to analyze
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
stack2Load  = 1;
path2Stacks = strcat(images2Analyze(stack2Load).folder,'\');
p2file      = strcat(path2Stacks,images2Analyze(stack2Load).name);
warning('off');
fileInfo    = loadMovie.tif.getinfo(p2file);
IM     = loadMovie.tif.getframes(p2file, 10); %Loading on of the frame
warning('on');

IM = imgaussfilt(IM,2);
%Adaptive threshold (usually better for high conc)
T = adaptthresh(IM,Threshold,'ForegroundPolarity','dark','Statistic','mean');
BWadapt = imbinarize(IM,T);
BWadapt = bwareaopen(BWadapt,4);
se = strel('disk',2);
BWadapt = imclose(BWadapt,se);

%No Adaptive Threshold (usually better for low conc)
BW = imbinarize(IM);
holes = ~BW;
holes = bwareaopen(holes,9);

IM2 = IM;
IM3 = IM;
holesAdapt = ~BWadapt;
% IM3(~BWadapt) = mean(mean(IM(~BWadapt))); 
% IM4 = IM2;
% IM4(IM4>0) = 1;

H0 = figure;
hold(gca,'on')
subplot(1,3,1)

imagesc(IM)
title('Raw image')
axis image

subplot(1,3,2)
imagesc(holesAdapt)
title('Adaptive Threshold')
axis image

subplot(1,3,3)
imagesc(holes)
title('Automated Threshold')
axis image
hold off

fileName0 = sprintf('%s%s%s-Pores_Threshold-%d',mainFolderName,'\',currentFolderName,Threshold*100);
savefig(H0,fileName0)

%TODO: allow User to change threshold (e.g. sliding bar?)
%% Asking for User input
prompt = {'Based on the Figure shown, do you want to use adaptive threshold (yes/no)'};
dlgTitle = 'Method to use for the analysis';
numLines = 1;
defaultVal = {'yes'};
answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);

%% testing and storing user input
assert(~isempty(answer{1}),'User cancel dialog, calculation aborted');
assert(ischar(answer{1}), 'Unexpected answer, valid answer are: yes, no, y, n');

if (or(strcmp(answer{1},'yes'),strcmp(answer{1},'y')))
    method2Use = 'adaptThreshold';
elseif or(strcmp(answer{1},'no'),strcmp(answer{1},'n'))
    method2Use = 'normThreshold';
else
    error('Unexpected answer, valid answer are: yes, no, y, n');
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
        I  = imgaussfilt(I,3);
        
        if strcmp(method2Use,'adaptThreshold')
        % Get an adaptive threshold (~depends on local intensity levels)
        T  = adaptthresh(I,Threshold,'ForegroundPolarity','dark','Statistic','mean');
        % Binarization of image and cleaning
        BWadapt = imbinarize(I,T);
        BWadapt = bwareaopen(BWadapt,4);
        se = strel('disk',2);
        BWadapt = imclose(BWadapt,se);
        BW_pores = ~BWadapt;
        
        elseif strcmp(method2Use,'normThreshold')
        BW = imbinarize(IM2proc);
        BW_pores = ~BW;
        BW_pores = bwareaopen(BW_pores,9); 
        end
        
        %Remove pores that are in touch with the edges(uncertainty about
        %size)
        %Need to change this so it does not always clear pores
        BW_pores = imclearborder(BW_pores);
        
        [L,n] = bwlabel(BW_pores);

        outData = zeros(n,3);
        for k = 1:n
        tmpBW = L==k;
        [ fWidth ] = SDcalc.fastWidthBW( tmpBW );
        tmpSize = sum(sum(tmpBW));
        outData(k,1) = fWidth;
        outData(k,2) = tmpSize;

        [B] = bwboundaries(tmpBW,'noholes');
        Blength = cellfun(@length,B);
        [~, idx] = maxk(Blength,2);
        B = B(idx);
        assert(length(B)==1,'problems');
        boundary = B{1};
        [ vals, names ] = SDcalc.solidity( boundary' );
        outData(k,3) = vals{1};

        end
        
%         %Get properties of the pores on the image.
%         stats = regionprops(BW_pores,'Area','MajorAxisLength',...
%         'MinorAxisLength');

%         dataForHistogram = [dataForHistogram [stats.Area].*pxArea];
    end
    
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

