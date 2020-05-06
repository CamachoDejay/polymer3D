function [frameInfo,movInfo] = getInfo(path2File)
    
    [path,file,ext] = fileparts(path2File);
    imReader = VideoReader(path2File);
    
    %get the number of frame
    maxFrame = round(imReader.Duration*imReader.FrameRate);
    %get the exposure time
    expT     = 1/imReader.FrameRate*1000; %in ms
    %store a few thing for the rest of the code
    isMultiImage = false;
    
    isZStack = false;
    
    Cam  = 0;
    
    X = imReader.Width;
    Y = imReader.Height;
        
    %Store info for output
    frameInfo.File = [file ext];
    
    movInfo.Width  = X;
    movInfo.Length = Y;
    movInfo.Path   = path;
    movInfo.isMultiImage = isMultiImage;
    movInfo.isZStack = isZStack;
    movInfo.Cam = Cam;
    movInfo.expT = expT;
    movInfo.maxFrame = maxFrame;