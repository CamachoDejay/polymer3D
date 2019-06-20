I function [ frameCell, movieInfo, tfl ] = getInfo( path2file )
%GETINFO receives as input the path to the file and gives back all
%infromation about the frames, the movie and total number of frames
%   Detailed explanation goes here
warning('off','all')
tObj = Tiff(path2file,'r');

movieInfo.Width  = tObj.getTag(256);
movieInfo.Length = tObj.getTag(257);
idx = strfind(path2file,filesep);
movieInfo.Path   = path2file(1:(idx(end))-1);

assert(tObj.currentDirectory == 1)
header = tObj.getTag(270);

% test here for having multiple image files
imFileStart = strfind(header, '<Image');
imFileEnd   = strfind(header, '</Image>');
nImFiles = length(imFileStart);
assert(nImFiles == length(imFileEnd), 'Problems with image files, size does not match')

if nImFiles > 1
    movieInfo.isMultiImage = true;
else
    movieInfo.isMultiImage = false;
end

frameCell = cell(nImFiles,1);

for imIdx = 1:nImFiles

    % get frame header
    frameHead = header(imFileStart(imIdx):imFileEnd(imIdx)+7);
    % Index frame header to find information relative to each tif index
    [k1, k2, k3, k4, nFrames] = indexFrameHeader(frameHead);
    % init frame info structure
    frameInfo = initFrameInfoStruc(nFrames);
    % fill in frame info
    for i = 1:nFrames
        str1 = frameHead(k1(i):k2(i)-1);
        str2 = frameHead(k3(i):k4(i));
        [ frameInfo(i).C, frameInfo(i).T, frameInfo(i).Z, frameInfo(i).IFD,...
            frameInfo(i).P, frameInfo(i).File, frameInfo(i).Pos,...
            frameInfo(i).expT,frameInfo(i).time] = getInfoFromString( str1, str2 );
 
    end
    tZero = frameInfo(1).time;
    for i = 1:nFrames
       frameInfo(i).time = frameInfo(i).time-tZero; 
    end
        
    frameCell{imIdx} = frameInfo;
end
% check frameInfo

[checkRes] = checkFrameInfo(frameInfo);
%if checkRes is true we fix camera frame
if checkRes
    %TODO: Fix misynchronization
    error('Fixing synchronization is not ready yet, sorry for the inconvenience')
end

%Add extrainfo to the movie, in particular, info about camera, max Frame,
%and zStack into movieInfo
movieInfo.isZStack = size(unique({frameInfo.T}),2)==1;
movieInfo.Cam      = str2double(unique({frameInfo.C}));
movieInfo.expT = frameInfo(1).expT;
switch movieInfo.isZStack
    case 0
        for i = 1: size(movieInfo.Cam,2)
            idx2Cam = strcmp({frameInfo(:).C},num2str(movieInfo.Cam(i)));
            maxFrame = max(str2double({frameInfo(idx2Cam).T}))+1;
            movieInfo.maxFrame(i) = maxFrame;
        end
    case 1
         for i = 1: size(movieInfo.Cam,2)
            idx2Cam = strcmp({frameInfo(:).C},num2str(movieInfo.Cam(i)));
            maxFrame = max(str2double({frameInfo(idx2Cam).Z}))+1;
            movieInfo.maxFrame(i) = maxFrame;
        end
end


tfl = 0; % Total frame length
while true
    tfl = tfl + 1; % Increase frame count
    if tObj.lastDirectory(), break; end
    tObj.nextDirectory();
end

tObj.setDirectory(1)
warning('on','all')

tObj.close

% this is so we dont break backwards compability, but we might want to
% change this into something cleaner later
if ~movieInfo.isMultiImage
    % then we should behave as before, thus we overwrite the cell by the
    % structure.
    frameCell = frameInfo;
    
end

end



function [k1, k2, k3, k4, nFrames] = indexFrameHeader(frameHeader)
% helper function that indexes over the frame header
    k1 = strfind(frameHeader, '<TiffData');
    k2 = strfind(frameHeader, '</TiffData>');
    k3 = strfind(frameHeader, '<Plane');
    k4 = strfind(frameHeader, '"/>');
    k4(k4<min(k3)) = [];
    nFrames = size(k3,2);

    if(~all([length(k2), length(k3), length(k4)]==nFrames))
        warning('Trying to find error and correct file (exposure time was likely too low or area imaged too large');
        if length(k1)>length(k3)
            idx = strfind(frameHeader,'IFD=');
            idx1 = strfind(frameHeader, 'FirstT=');
            idx2 = strfind(frameHeader, 'FirstC=');
            idx3 = strfind(frameHeader,'PlaneCount=');
            idx4 = strfind(frameHeader,'FirstZ=');
           
            IFD = zeros(size(idx1));
            C   = IFD;
            T   = IFD;
            
            for i = 1:size(idx1,2)
                IFD(i) = str2double(frameHeader(idx(i)+5:idx3(i)-3));
                C(i)   = str2double(frameHeader(idx2(i)+8:idx1(i)-3));
                T(i)   = str2double(frameHeader(idx1(i)+8:idx4(i)-3));
            end
       
            %find which frame is duplicated
            val = unique(T); % which will give you the unique elements of A in array B
            Ncount = histc(T, val);
            
            duplicate = val(Ncount>2);
            nWrong = length(duplicate);
            disp([num2str(nWrong) ' wrong lines found']);
            if nWrong > 0.01*nFrames
                error('More than 1% of the data has mistakes, cannot pursue')
            end
            
            idx2Delete = zeros(1,nWrong);
            %correct the errors
            for i = 1:nWrong
                currErr = duplicate(i);
                cIFDs = IFD(T==currErr);
                cCs   = C(T==currErr);
                
                testCam = sum(cCs);
                %determine which camera is wrong
                if testCam > 1
                    cIFDs = cIFDs(cCs==1);
                    cCam = 1;
                else
                    cIFDs = cIFDs(cCs==0);
                    cCam = 0;
                end
                
                %determine which of these IFD value is duplicated
                counter = zeros(size(cIFDs));
                for j = 1:length(cIFDs)
                    counter(j) = sum(IFD==cIFDs(j));
                     
                end
                
                [~,idx] = max(counter);
                
                val2Delete = cIFDs(idx);
                idx2Delete(i) = find(and(C==cCam,and(IFD==val2Delete,T==duplicate(i))));

                
            end
            k1(idx2Delete) = [];
            k2(idx2Delete)  = [];
            
            
        end
        
    end

end

function out = initFrameInfoStruc(nFrames)
% helper function to init an empty frameInfo structure
    out(nFrames).C    = [];
    out(nFrames).T    = [];
    out(nFrames).Z    = [];
    out(nFrames).IFD  = [];
    out(nFrames).P    = [];
    out(nFrames).File = [];
    out(nFrames).Pos  = [];
    out(nFrames).expT = [];
    out(nFrames).time  = [];
end

function [checkRes] = checkFrameInfo(frameInfo)
    disp('checking Camera synchronization');
    frame2Comp = 20;
    cellC = {frameInfo.C};
    matC = cellfun(@str2num,cellC);
    %We check the first 20 frames as they should be perfectly synchronized
    %if camera sync was properly used.
    test = abs(diff(matC(1:frame2Comp)));
    checkRes = false;
    if all(test)
        
    else
        answer = questdlg('It seems like the camera are not properly synchronized, do you still want to proceed?','Question to user','No','Yes','No');
        switch answer
            case 'Yes'
                checkRes = true;
            case 'No'
                disp('If you are running folder analysis, please remove the unsynced file from the folder');
                error('Camera are not synchronized, User aborted the analysis');
            otherwise
                disp('If you are running folder analysis, please remove the unsynced file from the folder');
                error('Camera are not synchronized, User aborted the analysis');
                
        end
    end
    
    

end