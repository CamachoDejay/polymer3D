function [frame] = getFrame(path,idx)
assert(isfile(path),'path given is not a file');
[~,~,ext] = fileparts(path);
assert(strcmpi(ext,'.his'),'.his file requested to load frame');

nFrames = length(idx);
%Create a bioformat reader
imReader = bfGetReader(path);

frame =zeros(imReader.getSizeY,imReader.getSizeX,nFrames);
%data is read as series of image that contain only one plane so we first
%need to set the series to the correct frame and then get the current
%plane:
for i = 1: length(idx)
    
    fi = idx(i);
    imReader.setSeries(fi-1); %idx start from 1 in matlab but from 0 in Java

    frame = bfGetPlane(imReader,1); % get the first and only plane in the series
    %which for some reason needs to be one not 0 x)


end