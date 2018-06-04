function [ frameInfo, movieInfo, tfl ] = getInfo( path2file )
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
Idesc = tObj.getTag(270);

k1 = strfind(Idesc, '<TiffData');
k2 = strfind(Idesc, '</TiffData>');

k3 = strfind(Idesc, '<Plane');
k4 = strfind(Idesc, '"/>');
k4(k4<min(k3)) = [];


il = size(k1,2);
assert(all([length(k2), length(k3), length(k4)]==il),'Plane and tif information do not match');

frameInfo(il).C = [];
frameInfo(il).T = [];
frameInfo(il).Z = [];
frameInfo(il).IFD = [];
frameInfo(il).P = [];
frameInfo(il).File = [];
frameInfo(il).Pos = [];


for i = 1:il
    str1 = Idesc(k1(i):k2(i)-1);
    str2 = Idesc(k3(i):k4(i));

    [ frameInfo(i).C, frameInfo(i).T, frameInfo(i).Z, frameInfo(i).IFD, frameInfo(i).P, frameInfo(i).File, frameInfo(i).Pos ] = getInfoFromString( str1, str2 );
end

%Add info about camera, max Frame, and zStack into movieInfo

movieInfo.isZStack = size(unique({frameInfo.T}),2)==1;
movieInfo.Cam      = str2double(unique({frameInfo.C}));

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
end

