function [ frameInfo, movieInfo, tfl ] = getInfo( path2file )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
warning('off','all')
tObj = Tiff(path2file,'r');

movieInfo.Width  = tObj.getTag(256);
movieInfo.Length = tObj.getTag(257);
idx = strfind(path2file,filesep);
movieInfo.Path   = path2file(1:(idx(end))-1);

assert(tObj.currentDirectory == 1)
Idesc = tObj.getTag(270);

k1 = strfind(Idesc, '<TiffData');
k2 = strfind(Idesc, '</TiffData>');

il = size(k1,2);
frameInfo(il).C = [];
frameInfo(il).T = [];
frameInfo(il).Z = [];
frameInfo(il).IFD = [];
frameInfo(il).P = [];
frameInfo(il).File = [];


for i = 1:il
    str = Idesc(k1(i):k2(i)-1);

    [ frameInfo(i).C, frameInfo(i).T, frameInfo(i).Z, frameInfo(i).IFD, frameInfo(i).P, frameInfo(i).File ] = getInfoFromString( str );
end

%Add info about camera, max Frame, and zStack into movieInfo

movieInfo.isZStack = size(unique({frameInfo.T}),2)==1;
movieInfo.Cam      = str2double(unique({frameInfo.C}));

switch movieInfo.isZStack
    case 0
        for i = 1: size(movieInfo.Cam,2)
            idx2Cam = strcmp({frameInfo(:).C},num2str(movieInfo.Cam(i)));
            maxFrame = max(str2double({frameInfo(idx2Cam).T}));
            movieInfo.maxFrame(i) = maxFrame;
        end
    case 1
         for i = 1: size(movieInfo.Cam,2)
            idx2Cam = strcmp({frameInfo(:).C},num2str(movieInfo.Cam(i)));
            maxFrame = max(str2double({frameInfo(idx2Cam).Z}));
            movieInfo.maxFrame(i) = maxFrame;
        end
end






tfl = 0; % Total frame length
while true
    tfl = tfl + 1; % Increase frame count
    if tObj.lastDirectory(), break; end;
    tObj.nextDirectory();
end

tObj.setDirectory(1)
warning('on','all')

tObj.close
end

