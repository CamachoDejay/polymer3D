function [ movC1, movC2, idx ] = loadMulticamOME( frameInfo, movieInfo, frames )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

assert(min(size(frames))==1,'Frames must be a vector')
assert(isstruct(frameInfo),'info must be a structure')
assert(max(frames)<= length(frameInfo) / 2, 'your request exceeds the total number of frames')
ImL = movieInfo.Length;
ImW = movieInfo.Width;
path2omes = movieInfo.Path;


nFramres = length(frames);
idx      = zeros(nFramres,2);
Cinfo    = {frameInfo.C};
Finfo    = {frameInfo.T};


for cam = 0:1
    C = num2str(cam);
    for fi = 1:nFramres
        frame = frames(fi);
        F = num2str(frame-1);
        tmp = and(strcmp(Cinfo,{C}),strcmp(Finfo,{F}));
        test = find(tmp);
%             if length(test)>1
%                 disp('HAAAAA')
%             elseif isempty(test)
%                 disp('NOOOO')
%             end
        idx(fi,cam+1) = find(tmp);
    end
end


movC1 = uint16(zeros(ImL,ImW,nFramres));
movC2 = uint16(zeros(ImL,ImW,nFramres));

h = waitbar(0,'Please wait...');
steps = nFramres*2;

oldPath = '';
warning('off','all')
for fi = 1:nFramres
    i      = idx(fi,1);
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
    movC1(:,:,fi) = tObj.read;
    waitbar(fi / steps)
end

for fi = 1:nFramres
    i      = idx(fi,2);
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
    movC2(:,:,fi) = tObj.read; 
    waitbar((fi + nFramres) / steps)
end
tObj.close
close(h) 
warning('on','all')

end

