clear;
clc;
close all;
%% User Input
pxSize = 100; %in nm
Threshold = 0.18; %number between 0 and 1 (% of max intensity)
bigPores = 50; %in px ("draw a line" between small pores and big pores);

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
Folders  = Folder_Content(3:end,:);
folderContent = Folders([Folders.isdir]==1);

for j = 1:size(folderContent,1)
    currentFolderPath = strcat(folderContent(j).folder,'\',folderContent(j).name);
    currentFolderContent = dir(currentFolderPath);
    
    index2Images      = contains({currentFolderContent.name},'.tif');
    images2Analyze = currentFolderContent(index2Images);
    path2Images = strcat(images2Analyze(1).folder,'\');
    
    if j == 1
        % Loading image 1
        p2file = strcat(path2Images,images2Analyze(1).name);
        fileInfo = loadMovie.tif.getinfo(p2file);
        IM = loadMovie.tif.getframes(p2file, 1);


        % Displaying the First image of stack
        figure(1)
        imagesc(IM(:,:,1))
        colormap('hot')

        % Image Processing, Histogram
        I = imgaussfilt(IM, 2);
        % double to do division
        I = double(I);
        I = I./max(max(I));

        vIm = I(:);
        figure(2)
        histogram(vIm);
        % Binary Image and Threshold

        %extract the network as a binary image
        BW_network = I>Threshold;
        %extract the pores
        BW_pores = I<Threshold;
        BW_pores = imclearborder(BW_pores);

        %plot the Figure
        figure(3)
        subplot(1,3,1)
        imagesc(I)
        colormap('hot')
        axis image

        subplot(1,3,2)
        imagesc(BW_network)
        % colormap('hot')
        axis image

        subplot(1,3,3)
        imagesc(BW_pores)
        axis image
    end
% Loop through all files

    Results = struct('numPores',cell(1,size(images2Analyze,1)),'numBigPores',cell(1,size(images2Analyze,1)),'Area',...
        cell(1,size(images2Analyze,1)),'areaBigPores',cell(1,size(images2Analyze,1)),...
        'meanArea',cell(1,size(images2Analyze,1)),'meanAreaBigPores',cell(1,size(images2Analyze,1)),...
        'medArea',cell(1,size(images2Analyze,1)),'medAreaBigPores',cell(1,size(images2Analyze,1)),...
        'MajorAxis',cell(1,size(images2Analyze,1)),'MinorAxis',cell(1,size(images2Analyze,1)),...
        'meanElip',cell(1,size(images2Analyze,1)),'medElip',cell(1,size(images2Analyze,1)));

    h = waitbar(0, 'Analyzing images');
    for i=1:size(images2Analyze,1)

        %loading image number i
        p2file = strcat(path2Images,images2Analyze(i).name);
        fileInfo = loadMovie.tif.getinfo(p2file);
        IM = loadMovie.tif.getframes(p2file, 1);
        I = imgaussfilt(IM, 2);
        I = double(I);
        I = I./max(max(I));

        %extract the network as a binary image
        BW_network = I>Threshold;
        %extract the pores
        BW_pores = I<Threshold;
        BW_pores = imclearborder(BW_pores);

        stats = regionprops(BW_pores,'Area','MajorAxisLength',...
        'MinorAxisLength');
        
        if isempty(stats) %avoid concatenation errors due to empty cell later
            clear stats;
            stats.Area = NaN;
            stats.MajorAxisLength = NaN;
            stats.MinorAxisLength = NaN;
            Results(i).numPores = NaN;
            
        else
        Results(i).numPores = size(stats,1);
        Results(i).Area = [stats.Area].*pxArea;
        Results(i).MajorAxis = [stats.MajorAxisLength].*pxArea;
        Results(i).MinorAxis = [stats.MinorAxisLength].*pxArea;
        end
  
        Results(i).numBigPores = size(stats([stats.Area]>bigPores),1);
        Results(i).areaBigPores = Results(i).Area(Results(i).Area>bigPores);
        
        if isempty(Results(i).areaBigPores)%avoid concatenation errors due to empty cell later
            Results(i).areaBigPores = 0;
        end
        
        Results(i).meanArea = mean(Results(i).Area);
        Results(i).meanAreaBigPores = mean(Results(i).Area(Results(i).Area>bigPores));
        Results(i).medArea  = median(Results(i).Area);
        Results(i).medAreaBigPores = median(Results(i).Area(Results(i).Area>bigPores));
       
        Results(i).Ellipticity = [stats.MinorAxisLength]./[stats.MajorAxisLength];
        Results(i).meanElip = mean(Results(i).Ellipticity);
        Results(i).medElip = median(Results(i).Ellipticity);
        waitbar(i/size(images2Analyze,1));
    end
    close(h);
    %% Save results to excel sheets

    fileNameExcel = sprintf('%s%s%s-ResultsSummary',mainFolderName,'\',...
        currentFolderName);
    matName =  regexprep(folderContent(j).name,'\.','_');
    fileNameMat   = sprintf('%s%s%s-FullResults',folderContent(j).folder,'\',...
        matName);
    save(fileNameMat,'Results');

    Title=[ 'numPores    ';'numBigPores ';'meanArea    ';'meanBigPores';...
        'medArea     ';'medBigPores '; 'meanEllip   ';'medEllip    '];
        Titletable=cellstr(Title)';
        xlswrite(fileNameExcel,Titletable,folderContent(j).name,'A1:H1');

    Result2print = [[Results.numPores]', [Results.numBigPores]',...
        [Results.meanArea]',[Results.meanAreaBigPores]',[Results.medArea]',...
        [Results.medAreaBigPores]',[Results.meanElip]',[Results.medElip]'];
        sizeRes = size(Results,2);
        Range = sprintf('A2:H%d',sizeRes+1);
        xlswrite(fileNameExcel,Result2print,folderContent(j).name,Range);

    
end
h = msgbox('The Data were succesfully saved !', 'Success');
