function [frame] = getFrame(path,idx)

    assert(isfile(path),'path given is not a file');
    [~,~,ext] = fileparts(path);
    assert(strcmpi(ext,'.mpg'),'.his file requested to load frame');

    nFrames = length(idx);
    
    imReader = VideoReader(path);
    
    
    frame = int32(zeros(imReader.Height,imReader.Width,nFrames));
    %data is read as series of image that contain only one plane so we first
    %need to set the series to the correct frame and then get the current
    %plane:
    for i = 1: length(idx)

        fi = idx(i);
        frame(:,:,i) = int32(rgb2gray(read(imReader,fi))); % get the first and only plane in the series
        %which for some reason needs to be one not 0 x)


    end