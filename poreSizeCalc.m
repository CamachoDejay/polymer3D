clear;
clc;
close all;
%% User Input
prompt = {'Enter the pixel size: ','Enter number of frame to analyze: '};
dlgTitle = 'User input for Pore size calculation';
numLines = 1;
defaultVal = {'100','244'};
answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);

%% Checking user input
assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')

pxSize = str2double(answer(1));
assert(~isnan(pxSize),'Number of Frame should be numerical');%If not a number

nFrame = str2double(answer(2));
assert(~isnan(nFrame),'Number of Frame should be numerical');%If not a number

fileExt = '.tif';
outputName = 'PoreSize-Results';

%% Loading Data
% conversion from pixel to area
pxArea = pxSize*pxSize*1e-6; %in µm^2

[file2Analyze,currentFolderName,outDir] = Load.Folder(fileExt,outputName);

%% Looping through the Data

h = waitbar(0);
nImStacks = size(file2Analyze,1);

allDataAdapt = [];
allDataAuto = [];
for j = 1:nImStacks
    hMessage = sprintf('Loading image stack number %d/%d',j,nImStacks);
    waitbar(0,h,hMessage);
    %Data loading
    path2Stacks = strcat(file2Analyze(j).folder,filesep);
    tmpName = file2Analyze(j).name;
    p2file      = strcat(path2Stacks,tmpName);
    warning('off','all')
    fileInfo    = Load.Movie.tif.getinfo(p2file);
    
    warning('on','all')
    tNframes = fileInfo.Frame_n;
    assert(tNframes>=nFrame,'Requested number of frame is larger than the number of frame in the file')
    
    % init data that contains all infor for a single tif file
    tifStackData = [];
    
    hMessage = sprintf('Analysis of Stack Number %d/%d',j,nImStacks);
    %loop through the frames of the current stack
    nIM = nFrame;
    waitbar(j/nImStacks,h,hMessage);
    
    for i=1:nIM %% PARFOR CAN BE PLACED HERE
        
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
        %TODO:
        %Calculate pore volume fraction
    end
    
    % add info about tif index
    tTmp = table(ones(size(tifStackData,1),1).*j,'VariableNames',{'TifIDX'});
    tifStackData = [tTmp, tifStackData];

    disp('---------------------NEXT TIF ----------')
    
    % store in main table
    if ~isempty(strfind(file2Analyze(j).name,'adapt'))
        allDataAdapt = [allDataAdapt; tifStackData];
    else
        allDataAuto = [allDataAuto; tifStackData];
    end
end
close(h);

%% saving data
% % % T = array2table(allData,...
% % %     'VariableNames',{'TifIDX','ImageIDX','Width','Area', 'Solidity','IsAtBorder'});

infoFileName = [outDir filesep 'info.txt'];
    fid = fopen(infoFileName,'wt');
    fprintf(fid,'This file contains information intended for the user of poreSizeCalc\n');
    fprintf(fid,' In such a way that the user knows what variable value were used.\n\n');
    fprintf(fid,'Pixel size used: %d',pxSize);
    fprintf(fid,'Number of frame analyzed: %d/%d\n',nFrame,tNframes);
    fclose(fid);

save([ outDir 'Adaptive-poreProps.mat'],'allDataAdapt')
save([ outDir 'Automated-poreProps.mat'],'allDataAuto')
h = msgbox('The Data were succesfully saved !', 'Success');

%% Plotting
totalAdapt = allDataAdapt.Area;
totalAdapt = totalAdapt(:);
[CDF,CCDF] = Plotting.getCDF(totalAdapt);

totalGlobal = allDataAuto.Area;
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
