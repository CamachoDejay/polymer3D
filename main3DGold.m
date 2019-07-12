%% GOLD TRACKiNG
clear 
clc
close all

%% User input

path2Cal  = 'F:\Data\Leuven Data\2019\07 - July\20190708\2D';
file.path = 'F:\Data\Leuven Data\2019\07 - July\20190708';
file.ext  = '.ome.tif'; 
delta = 50;% in px Size of the ROI around particles detected(radius 50 = 100x100 pixel
nParticles = 2;%number of particles expected in the movie has to be exact
pxSize = 95;%in nm
minDist = 6; %in pixels (min distant expected between particles
scaleBar = 2; %in um
tail = 20;%Length of the tail in frames, for plotting the traces on top of the movie
frameRate = 30; %for saving the movie
info.type = 'Transmission';%Transmission or normal movie
info.runMethod = 'load';
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
mov.g2.showFrame(1,5);

%% detection of the center of the beads
%preallocate memory for storing data of the different files
allData = struct('fileName',[],'locPos',[]);
allData(size(folder2Mov,2)).locPos = [];

for i =1: size(folder2Mov,2)
    
    currentPath = folder2Mov(i).folder;
    switch toAnalyze
        case '.mp4'
            p2file = [folder2Mov(i).folder filesep folder2Mov(i).name];
            v = VideoReader(p2file);%Create VideoReader object
            nFrames = floor(v.Duration*v.FrameRate);%extract the number of frames
            fullStackIn = zeros(v.Height,v.Width,nFrames);%preallocate memory
            for j = 1:nFrames
                frame = readFrame(v);%read frames

                fullStackIn(:,:,j) = rgb2gray(frame);%extract the intensity data from frames
            end
        otherwise
            myMov = Core.Movie(currentPath,info,'.tif');%Create Movie Object
            fullStack = myMov.getFrame(1);%extract first frame
            frame = fullStack.Cam1;
            %check if cropping is necessary
            if size(frame,2) > 400
                myMov.cropIm;
            end
            %load full stack
            fullStack = myMov.getFrame;
            %revert the intensity scale
            if strcmpi(info.type,'transmission')
                fullStackIn = imcomplement(fullStack.Cam1);
            else
                fullStackIn = fullStack.Cam1;
            end
    end
    
    frame = 5;
    %get the n maxima where n is the number of particles expected and
    %minDist is the distance expected between them
    [pos] = goldProj.nMaxDetection (fullStackIn(:,:,frame),nParticles,minDist);
    
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



end