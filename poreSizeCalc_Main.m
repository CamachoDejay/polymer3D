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
index2Images   = contains({Folder_Content.name},'.tif');
images2Analyze = Folder_Content(index2Images);

h = waitbar(0);
for j = 1:size(images2Analyze,1)
    hMessage = sprintf('Loading image stack number %d/%d',j,size(images2Analyze,1));
    waitbar(0,h,hMessage);
    path2Stacks = strcat(images2Analyze(j).folder,'\');
    p2file      = strcat(path2Stacks,images2Analyze(j).name);
    warning('off');
    fileInfo    = loadMovie.tif.getinfo(p2file);
    frames2Load = 1:1:fileInfo.Frame_n;
    IMStack     = loadMovie.tif.getframes(p2file, frames2Load);
    warning('on');
    Results = struct('numPores',cell(1,size(images2Analyze,1)),'numBigPores',cell(1,size(images2Analyze,1)),'Area',...
        cell(1,size(images2Analyze,1)),'areaBigPores',cell(1,size(images2Analyze,1)),...
        'meanArea',cell(1,size(images2Analyze,1)),'meanAreaBigPores',cell(1,size(images2Analyze,1)),...
        'medArea',cell(1,size(images2Analyze,1)),'medAreaBigPores',cell(1,size(images2Analyze,1)),...
        'MajorAxis',cell(1,size(images2Analyze,1)),'MinorAxis',cell(1,size(images2Analyze,1)),...
        'meanElip',cell(1,size(images2Analyze,1)),'medElip',cell(1,size(images2Analyze,1)));

    hMessage = sprintf('Analysis of Stack Number %d/%d',j,size(images2Analyze,1));
    for i=1:size(IMStack,3)
    waitbar(i/size(IMStack,3),h,hMessage);
        %loading image number i
        IM = IMStack(:,:,i);
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
            Results(i).numPores = 0;
            
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
        
    end
    close(h);
    %% Save results to excel sheets

    fileNameExcel = sprintf('%s%sResultsSummary',mainFolderName,'\');
    matName =  regexprep(Folder_Content(j).name,'\.','_');
    fileNameMat   = sprintf('%s%s%s-FullResults',Folder_Content(j).folder,'\',...
        matName);
    save(fileNameMat,'Results');
    sheetName = sprintf('Stack%d',j);
    Title=[ 'numPores    ';'numBigPores ';'meanArea    ';'meanBigPores';...
        'medArea     ';'medBigPores '; 'meanEllip   ';'medEllip    '];
        Titletable=cellstr(Title)';
        xlswrite(fileNameExcel,Titletable,sheetName,'A1:H1');

    Result2print = [[Results.numPores]', [Results.numBigPores]',...
        [Results.meanArea]',[Results.meanAreaBigPores]',[Results.medArea]',...
        [Results.medAreaBigPores]',[Results.meanElip]',[Results.medElip]'];
        sizeRes = size(Results,2);
        Range = sprintf('A2:H%d',sizeRes+1);
        
        xlswrite(fileNameExcel,Result2print, sheetName,Range);
end
h = msgbox('The Data were succesfully saved !', 'Success');
