function [ frameInfo, movieInfo, tfl ] = getomeinfo( path2file )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
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


for i = 1:il;
    str = Idesc(k1(i):k2(i)-1);

    [ frameInfo(i).C, frameInfo(i).T, frameInfo(i).Z, frameInfo(i).IFD, frameInfo(i).P, frameInfo(i).File ] = MovieLoad.ome.getInfoFromOMEstr( str );
end
tic
warning('off','all')
tfl = 0; % Total frame length
while true
    tfl = tfl + 1; % Increase frame count
    if tObj.lastDirectory(), break; end;
    tObj.nextDirectory();
end
toc
tObj.setDirectory(1)
warning('on','all')

tObj.close
end

