function [ frameCell, movieInfo, tfl ] = getInfo( path2file )
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
%checkRes ='Yes';
switch checkRes
    case 'Yes'
    case 'Fix'
        %try
        [frameInfo] = fixCamTiming(frameInfo); 
        warning('on')
        warning('Found some synchronization issue and fixed them');
        warning('off')
%         catch
%             warning('on')
%             warning('Something went wrong when trying to fix camera sync, no fix was apply');
%             warning('off');
%         end
      %  error('Fixing synchronization is not ready yet, sorry for the inconvenience');
    case 'No'
        disp('If you are running folder analysis, please remove the file from the folder');
        disp(['The unsynced file is: ',movieInfo.Path]);
        error('Camera are not synchronized, user aborted the analysis')
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
            if nWrong > ceil(0.02*nFrames)
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
    frame2Comp = 10;
    if length(frameInfo) < frame2Comp
       
        frame2Comp = size(frameInfo,1);
        
    end
    cellC = {frameInfo.C};
    matC = cellfun(@str2num,cellC);
    %We check the first 20 frames as they should be perfectly synchronized
    %if camera sync was properly used.
    test = abs(diff(matC(1:frame2Comp)));
    
    if all(test)
        checkRes = 'Yes';
    else
        checkRes = 'Fix';
    end   
end

function [frameInfo] = fixCamTiming(frameInfo)
    tmpInfo = frameInfo;
    timing = [tmpInfo.time];
    %Check that the timing is correct by taking the difference. Where the
    %camera are synchronize it will give 0 and expTime if camera are
    %asynchronize, it will be exptime everywhere in the bingging and the
    %end
    test = diff(timing);
    %extract the indices of the desync part
    idxStart = find(test<mean(test),1,'first');
    idxEnd  = find(test<mean(test),1,'last');
    %delete the data
    tmpInfo([1:idxStart-1,idxEnd:end])=[];
    
    %check that we kept the same number of frame for the two cameras
    camera = cellfun(@str2num,{tmpInfo.C});  
    assert(sum(camera==0)==sum(camera==1),'Something went wrong when fixing the sync');
    
    newTiming = [tmpInfo.time];
    %let us renumber the T after the deletion
   if str2double(tmpInfo(1).T) ==0
        refCam= camera(1);
        modifier = 0;

    elseif str2double(tmpInfo(2).T) ==0
        refCam= camera(2);
        modifier = 1;
   end

    idxRefCam = find(strcmp({tmpInfo.C},num2str(refCam))==1);
    for i = 1: length(idxRefCam)
        currRefId = idxRefCam(i);
        currT = tmpInfo(currRefId).T;
        %get the time of the ref camera for the current frame
        timeCurrRef = tmpInfo(currRefId).time;

        %find another frame with similar timing
        id = find((abs([tmpInfo.time]-timeCurrRef))<4);
        assert(length(id) == 2, 'couldnt find 2 frames for that time point');

        id = id(id~= currRefId);
        tmpInfo(id).T = currT;
        
        %previous version:
%         %if current index is not reference camera we need to fix things
%         if camera(i) ~=refCam
%         
%             % prev: abs(newTiming(i)-newTiming)<mean(test)
%             tmpInfo(i).T = currT;
%             
%         else
%             %if it is the ref camera then currT becomes the camera
%             currT = tmpInfo(i).T;
%             
%         end
        
    end
    
    % test that the camera are indeed synchroneous
    maxT = str2double(tmpInfo(end).T);
    for i = 0:maxT
        
        idx = strcmp({tmpInfo.T},{num2str(i)});
        
        camDiff = {tmpInfo(idx).C};
        
        %test cam diff
        %#1 we test that there is only 2 camera frames
        assert(length(camDiff)<=2,'More than two camera for a single frame')
        %#2 
        assert(~strcmp(camDiff{1}, camDiff{2}),'The two cameras corresponding to the same time point are not different')
        
        %test timing difference
        camTiming = {tmpInfo(idx).time};
        
        assert(abs(camTiming{1} - camTiming{2}) < 4, 'Camera delay bigger than 4ms after fixing, something is wrong')
        
        
    end
    
    
    
    
    frameInfo = tmpInfo;
    warning('There was some issues with the synchronization, we had to fixed the out of sync frames');

end
