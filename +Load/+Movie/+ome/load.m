function [ movC1, movC2, idx ] = load( frameInfo, movieInfo, frames )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

assert(isvector(frames),'Frames must be a vector');
assert(isstruct(frameInfo),'info must be a structure');
assert(max(frames)<= min(movieInfo.maxFrame), 'your request exceeds the total number of frames');

isZStack = movieInfo.isZStack; % define if zStack or not
cam      = movieInfo.Cam; %give indexes of camera, size is number of cam
assert(~isempty(cam),'Camera information is missing.');
assert(and(isvector(cam),length(cam)<3),'Camera must be a vector of max. length 2.');
nFramres = length(frames);
Cinfo    = {frameInfo.C};

switch isZStack
    case 0
        Finfo = {frameInfo.T};
    case 1
        Finfo = {frameInfo.Z};
end

switch size(cam,2)
    case 1
        idx = zeros(nFramres,1);
    case 2 
        idx = zeros(nFramres,2);
    otherwise
        error('Unexpected number of cam loaded');
end
        
for iCam = 1:size(cam,2)
    C = num2str(cam(iCam));
    for fi = 1:nFramres
        frame = frames(fi);
        F = num2str(frame-1);
        tmp = and(strcmp(Cinfo,{C}),strcmp(Finfo,{F}));
        if all(tmp==0)
            disp('error');
        end
        idx(fi,iCam) = find(tmp);
    end
end

% 
warning('off','all')

assert (size(idx,2) == length(movieInfo.Cam), 'Problem with number of cameras');
%We always assume to have one camera
camIdx = idx(:,1);
[movC1,~] = loadSingleCam(movieInfo,frameInfo,frames,camIdx);

switch size(idx,2)
    
    case 1 %Only one Cam
         movC2 = []; %We give the same movie for both camera (or leave
        %empty?)
    case 2 %when two cameras are used
        camIdx = idx(:,2);
        [movC2,~] = loadSingleCam(movieInfo,frameInfo,frames,camIdx);  
        
    otherwise
        error('Your file has more cameras than we can handle atm (max. 2)');
        
end
warning('on','all')

end
function [mov,oldPath] = loadSingleCam(movieInfo,frameInfo,frames,camIdx)
    ImL      = movieInfo.Length; % Length of the frame
    ImW      = movieInfo.Width; %width of the Frame
    nFramres = length(frames); %frame is frame2load
    mov = uint16(zeros(ImL,ImW,nFramres));
    oldPath = '';
    path2omes = movieInfo.Path;
    h = waitbar(0,'Loading a camera');
    steps = nFramres;
    for fi = 1:nFramres
        i      = camIdx(fi);
        f2load = frameInfo(i).File;
        p2file = [path2omes filesep f2load];

        if ~strcmp(oldPath,p2file)

            if exist( 'tObj', 'var' )
                tObj.close
            end
            tObj   = Tiff(p2file,'r');
            oldPath = p2file;

        end

        TifDir = str2double( frameInfo(i).IFD ) + 1;
        %     disp(TifDir)
        tObj.setDirectory(TifDir)
        mov(:,:,fi) = tObj.read;
        waitbar(fi / steps)
    end
    close (h);
    tObj.close;
end
