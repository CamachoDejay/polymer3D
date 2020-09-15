function [ frameInfo, movInfo ] = getInfo( path2File )
%GETINFO Summary of this function goes here
%   Detailed explanation goes here
[path,file,ext] = fileparts(path2File);
warning('off','all')
tObj = Tiff(path2File,'r');

movInfo.Width  = tObj.getTag(256);
movInfo.Length = tObj.getTag(257);
movInfo.Path   = fileparts(path2File);

tfl = 0; % Total frame length
while true
    tfl = tfl + 1; % Increase frame count
    if tObj.lastDirectory(), break; end;
    tObj.nextDirectory();
end
tObj.setDirectory(1)
warning('on','all')

%store a few thing for the rest of the code
isMultiImage = false;

isZStack = false;
Cam  = 0;
%get the number of frame
maxFrame = tfl;
%get the exposure time
expT     = NaN; %in ms

movInfo.isMultiImage = isMultiImage;
movInfo.isZStack = isZStack;
movInfo.Cam = Cam;
movInfo.expT = expT;
movInfo.maxFrame = maxFrame;



tObj.close
%Store info for output
frameInfo.File = [file ext];
end

