clear;
clc;
close all;
%% User Input
prompt = {'Enter the pixel size: ','Enter number of frame to analyze: ',...
    'number of frame to ignore beginning:', 'number of frame to ignore end:','Stack to ignore:'};
dlgTitle = 'User input for Pore size calculation';
numLines = 1;
defaultVal = {'100','244','0','0','0'};
answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);

%% Checking user input
assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')

pxSize = str2double(answer(1));
assert(~isnan(pxSize),'Number of Frame should be numerical');%If not a number

nFrame = str2double(answer(2));
assert(~isnan(nFrame),'Number of Frame should be numerical');%If not a number

bIgnore = str2double(answer(3));
assert(~isnan(bIgnore),'Number of Frame to ignore should be numerical');

aIgnore = str2double(answer(4));
assert(~isnan(aIgnore),'Number of Frame to ignore should be numerical');

sIgnore = str2double(answer(3));
assert(~isnan(bIgnore),'Number of Frame should be numerical');

fileExt = '.tif';
outputName = 'PoreSize-Results';

%% Loading Data
% conversion from pixel to area
pxArea = pxSize*pxSize*1e-6; %in µm^2

[file2Analyze,currentFolderName,outDir] = Load.Folder(fileExt,outputName);

%% Looping through the Data

%h = waitbar(0);
nImStacks = size(file2Analyze,1);
idx2Stack = 1:nImStacks;
if sIgnore~=0
    nImStacks = nImStacks-length(sIgnore);
    idx2Stack(idx2Stack == sIgnore) = [];
end

if bIgnore ~=0
    startIdx = bIgnore;
else
    startIdx = 1;
end

if aIgnore ~= 0
    endIdx = nFrame-aIgnore;
else
    endIdx = nFrame;  
end

allDataAdapt = struct('filename', [], 'Data', []);
allDataAuto = struct('filename', [], 'Data', []);

allDataAdapt(length(idx2Stack)/2).filename = [];
allDataAuto(length(idx2Stack)/2).filename = [];
nAdapt = 1;
nAuto  = 1;
for j = 1:length(idx2Stack)
    jdx = idx2Stack(j);
   % hMessage = sprintf('Loading image stack number %d/%d',j,nImStacks);
    %waitbar(0,h,hMessage);
    %Data loading
    path2Stacks = strcat(file2Analyze(jdx).folder,filesep);
    tmpName = file2Analyze(jdx).name;
    p2file      = strcat(path2Stacks,tmpName);
    warning('off','all')
    fileInfo    = Load.Movie.tif.getinfo(p2file);
    
    warning('on','all')
    tNframes = fileInfo.Frame_n;
    assert(tNframes>=nFrame,'Requested number of frame is larger than the number of frame in the file')
    
    % init data that contains all infor for a single tif file
    tifStackData = [];
    
   % hMessage = sprintf('Analysis of Stack Number %d/%d',j,nImStacks);
    %loop through the frames of the current stack
    nIM = nFrame;
   % waitbar(j/nImStacks,h,hMessage);
    disp('Loading Data');
    totVol = 0;
    filledVolume = 0;
    parfor i=startIdx:endIdx %% PARFOR CAN BE PLACED HERE
        
        % Loading image number i
        warning('off','all')
        IM     = Load.Movie.tif.getframes(p2file, i);
        warning('on','all')
        disp(['loaded frame: ' num2str(i)])
        
         % pores not conected to the border
        im_clean = imclearborder(IM);
        % pores conected to the border
        im_border = and(IM,~im_clean);
        
        % quick calculation of shape descriptors
        stats1 = regionprops('table',im_clean,'Area','MajorAxisLength','MinorAxisLength','Orientation','Solidity','Eccentricity');
        stats1.IsBorder = zeros(size(stats1,1),1);
        
        %TODO: Store the outside pore somehow to keep track of them ?
        stats2 = regionprops('table',im_border,'Area','MajorAxisLength','MinorAxisLength','Orientation','Solidity','Eccentricity');
        stats2.IsBorder = ones(size(stats2,1),1);

        regData = [stats1;stats2];
        
        disp(['Done for Image ' num2str(i) '/' num2str(nIM)])
        disp('---------------------NEXT FRAME ----------')

        % changing from pixel^2 to micron^2
        regData.Area = regData.Area.*pxArea;
        % changing form pixel to micron
        regData.MajorAxisLength = regData.MajorAxisLength.*(pxSize*1e-3);
        regData.MinorAxisLength = regData.MinorAxisLength.*(pxSize*1e-3);
        % add frame index
        tTmp = table(ones(size(regData,1),1).*i,'VariableNames',{'FrameIDX'});
        regData = [tTmp, regData];
        % store
        tifStackData = [tifStackData; regData];
        totVol = totVol + numel(IM);
        filledVolume  = filledVolume + sum(sum(~IM));
       
    end
    
    % add info about tif index
    tTmp = table(ones(size(tifStackData,1),1).*j,'VariableNames',{'TifIDX'});
    tifStackData = [tTmp, tifStackData];

    disp('Storing Data')
    
    % store in main table
    if ~isempty(strfind(file2Analyze(j).name,'adapt'))
        
        %allDataAdapt = [allDataAdapt; tifStackData];
        allDataAdapt(nAdapt).filename = file2Analyze(j).name;
        allDataAdapt(nAdapt).Data = tifStackData;
        allDataAdapt(nAdapt).Data.filledVolume(1) = filledVolume;%stillPX
        allDataAdapt(nAdapt).Data.totVolume(1)    = totVol;%stillPX
        nAdapt = nAdapt+1;
    else
        allDataAuto(nAuto).filename = file2Analyze(j).name;
        allDataAuto(nAuto).Data = tifStackData;
        nAuto = nAuto+1;
    end
    disp('---------------------NEXT TIF ----------')
end
%close(h);

%% saving data
% % % T = array2table(allData,...
% % %     'VariableNames',{'TifIDX','ImageIDX','Width','Area', 'Solidity','IsAtBorder'});

infoFileName = [outDir filesep 'info.txt'];
    fid = fopen(infoFileName,'wt');
    fprintf(fid,'This file contains information intended for the user of poreSizeCalc\n');
    fprintf(fid,' In such a way that the user knows what variable value were used.\n\n');
    fprintf(fid,'Pixel size used: %d',pxSize);
    fprintf(fid,'Number of frame analyzed: %d\n',nFrame);
    fclose(fid);

save([ outDir 'Adaptive-poreProps.mat'],'allDataAdapt')
save([ outDir 'Automated-poreProps.mat'],'allDataAuto')
h = msgbox('The Data were succesfully saved !', 'Success');

%% Plotting

totalAdapt = [];
totalGlobal = [];
for i = 1 : length(allDataAdapt)
    
    totalAdapt  = [totalAdapt ; allDataAdapt(i).Data.Area];
    totalGlobal = [totalGlobal ;  allDataAuto(i).Data.Area];

end

totalAdapt = totalAdapt(:);
[CDF,CCDF] = Plotting.getCDF(totalAdapt);


totalGlobal  = totalGlobal (:);
[CDF2,CCDF2] = Plotting.getCDF(totalGlobal);

figure()
subplot(1,2,1)
plot(CCDF.x,CCDF.y)
a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF for AREA',' for all tif files in folder using adaptive threshold'})
xlabel('Pore size')
ylabel('CCDF - prob [0 1]')
xlim([0 10^5])
a.FontSize = 14;

subplot(1,2,2)
plot(CCDF2.x,CCDF2.y)
a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF for AREA',' for all tif files in folder using automated threshold'})
xlim([0 10^5])
xlabel('Pore size')
ylabel('CCDF - prob [0 1]')
a.FontSize = 14;
