%% GOLD TRACKiNG
clear 
clc
close all

%% User input

path2Cal  = 'F:\Data\Leuven Data\2019\07 - July\2D-DarkField\Small';
file.path = 'F:\Data\Leuven Data\2019\07 - July\DarkfieldFewParticles\2P';
file.ext  = '.ome.tif'; 
focusPlane = 5;%=2 af
delta = 50;% in px Size of the ROI around particles detected(radius 50 = 100x100 pixel
nParticles = 2;%number of particles expected in the movie has to be exact
pxSize = 95;%in nm
minDist = 4; %in pixels (min distant expected between particles
scaleBar = 2; %in um
tail = 20;%Length of the tail in frames, for plotting the traces on top of the movie
frameRate = 30; %for saving the movie
info.type = 'normal';%Transmission or normal movie
info.runMethod = 'load';
info.calibrate = false;
toAnalyze = '.ome.tif';%accepted: .mp4, .ome.tif, folder. (folder that contain folders of .ome.tif.
outputFolder = 'Results';%name of the folder to output the results

%% Loading
%we get the zCalibration directory
folder2Mov = dir(file.path);
folder2Mov = folder2Mov(cell2mat({folder2Mov.isdir}));
%loop through the content of the directory
mov = struct();
count = 0;
for i = 3:size(folder2Mov,1)
    %Check if the directory
    folderPath = [folder2Mov(i).folder filesep folder2Mov(i).name];
    file2Analyze = Core.Movie.getFileInPath(folderPath,file.ext);

    if ~isempty(file2Analyze)
        count = count+1;
        filetmp.path = file2Analyze.folder;
        filetmp.ext  = file.ext;
        tmp = Core.MPMovie(filetmp , path2Cal,info);
        
        if count == 1
            tmp.giveInfo;
        else
            tmp.info = obj.trackMovies.(['mov' num2str(1)]).getInfo; 
        end
        tmp.calibrate;
        mov.(['g' num2str(i-2)]) = tmp;


    else

        warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);

    end

end
disp('=======> DONE ! <========')
mov.g1.showFrame(1,5);

%% detection of the center of the beads

%preallocate memory for storing data of the different files
allData = struct('fileName',[],'traces',[]);
allData(length(mov)).locPos = [];
fields = fieldnames(mov);
frame = 35;

for i = 1:length(mov)
    currMov = mov.(fields{i});
    %check detection
    fullStack = currMov.getFrame(frame);
    %revert the intensity scale
    if strcmpi(info.type,'transmission')
        fullStackIn = imcomplement(fullStack);
    else
        fullStackIn = fullStack;
    end

    frame = 5;
    %get the n maxima where n is the number of particles expected and
    %minDist is the distance expected between them
    [pos] = goldProj.nMaxDetection (fullStackIn(:,:,focusPlane),nParticles,minDist);

    x0 = pos(:,2);
    y0 = pos(:,1);

    pos = [y0 x0];
    %calculate center of mass of particles position for cropping later
    cropPos = round(mean(pos,1));
    %plot to check that the max detection worked
    figure
    imagesc(fullStackIn(:,:,frame))
    hold on
    plot(pos(:,2),pos(:,1),'r+')
    
    
    %% Fitting
    nFrames = currMov.raw.movInfo.maxFrame(1);
    %get the X Y dommain
    x = 1:size(fullStackIn,2);
    y = 1:size(fullStackIn,1);
    [domX,domY] = meshgrid(x,y);
    dom(:,:,1) = domX;
    dom(:,:,2) = domY;
    %preallocate memory
       
    data2Store = zeros(nFrames,3,nParticles);
    fitMov = zeros(100,100,nFrames);
    h = waitbar(0,'Fitting Data');%create waiting bar
    %unSortedData =[];
    planePos = currMov.calibrated.oRelZPos;
    for j = 1:nFrames
        currentFrame = double(currMov.getFrame(j));
        %inital detection of particles on currentFrame
        [pos] = goldProj.nMaxDetection (currentFrame(:,:,focusPlane),nParticles,minDist);
 
        x0 = pos(:,2);
        y0 = pos(:,1);
        %Multiple gaussian fitting occurs here
        [gPar,resnorm,res] = Localization.Gauss.MultipleFitting(currentFrame(:,:,focusPlane),x0,y0,dom,nParticles); 
        g = reshape(gPar(5:end),[2,nParticles]);
        g = g';
        
        nPart = size(x0,1);
        z = x0;
        %Get Z position
        for k = 1:nPart
            partData.col = x0(k);
            partData.row = y0(k);
            coord = round(g);
%            [Mag] = Core.MPLocMovie.getZPhasorMag(partData,minDist,currentFrame);
            ROI = currentFrame(coord(k,2)-minDist:coord(k,2)+minDist,...
                coord(k,1)-minDist:coord(k,1)+minDist,:);
                
            toFit = squeeze(sum(sum(ROI)));
            domain = 1:length(toFit);
            guess.sig = minDist/2;
            guess.mu  = focusPlane;
            [Res,fit] = SimpleFitting.gauss1D(toFit,domain,guess);

            z(k) = Res(2);
            if or(z(k)<min(domain),z(k)>max(domain))
                z(k)   = NaN;                           
            else
            tmpZ = floor(z(k));
            fracZ = z(k)-tmpZ;
            z(k) = planePos(tmpZ)+fracZ*(planePos(tmpZ+1) - planePos(tmpZ));
            z(k) = z(k)*1000;
            end
        end
        
        
       % unSortedData = [unSortedData; g z ones(size(g,1),1)*j  ];
        
        %Generate an image of the Fit to be able to plot in case we want to
        %check
        xHRes = linspace(1,size(fullStackIn,2),100);
        yHRes = linspace(1,size(fullStackIn,1),100);
        [domHResX,domHResY] = meshgrid(xHRes,yHRes);
        domHRes(:,:,1) = domHResX;
        domHRes(:,:,2) = domHResY;
        F = Localization.Gauss.MultipleGauss(gPar, domHRes,nParticles);
        
        if j>1
            %Tracking based on MSD minimization 
            newOrder = goldProj.simpleTracking(gPar(5:end),prevPos);
            %reshaping to format the final data and sorting with new order
            gPos = reshape(gPar(5:end),[2,nParticles]);
            gPos = gPos(:,newOrder);
            prevPos = reshape(gPos,[1,nParticles*2]);
            gPos = reshape(gPos,[1,2,nParticles])*pxSize;
            gPos(:,3,:) = z;
            %store data
            data2Store(j,:,:) = gPos;
            clear gPos

        else
            %First frame we just reshape
            gPos = reshape(gPar(5:end),[2,nParticles]);
            gPos = reshape(gPos,[1,2,nParticles]);
            gPos(:,3,:) = z;
            data2Store(j,:,:) = gPos;
            %store the gPar in the prev for checking particle order
            prevPos = gPar(5:end);
            clear gPos
        end
        %store fittings
        fitMov(:,:,j) = F;
        %update waitbar value
        waitbar(j/nFrames,h,'Fitting Data')
    end
    
    %clear waitbar
    close(h);
    
    currentPath = currMov.raw.movInfo.Path;
    %save data to the current folder being analyze
    filename = [currentPath filesep 'LocalizationData.mat'];
    save(filename,'data2Store');
    %store data in allData
    allData(i).traces = data2Store;
    allData(i).fileName = currentPath;  
    %store fittings
    fitMov(:,:,j) = F;

end

%% convert Data to table


for i =1: size(allData,2)
    traces = cell(nParticles,1);
    
    for j = 1: nParticles
        tabData = table(zeros(nFrames,1),zeros(nFrames,1),zeros(nFrames,1),...
            zeros(nFrames,1),'VariableNames',{'row','col','z','t'});
        
        tabData.col = allData(i).traces(:,1,j);
        tabData.row = allData(i).traces(:,2,j);
        tabData.z = allData(i).traces(:,3,j);
        tabData.t = (1:nFrames)';
        
        traces{j} = tabData;
    end
    allData(i).traces = traces;
end
trackRes = allData;
filename = [file.path filesep 'trackRes.mat'];
save(filename,'trackRes');
h = msgbox('Data succesfully saved');