function [frameInfo,movInfo] = getInfo(path2File)
    %TODO: Extract exposure time 
    [path,file,ext] = fileparts(path2File);
    fid = fopen(path2File,'r'); 
    

    header = fread(fid,2050,'uint16=>uint16'); % 2050 uint16 = 4100 bytes = 32800 bits

    xDim = double(header(22));
    yDim = double(header(329));
    nFrames = double(header(724));
      
    %get the number of frame
    maxFrame = nFrames;
    %get the exposure time
    expT     = NaN; %in ms
    %store a few thing for the rest of the code
    isMultiImage = false;
    
    isZStack = false;
    
    Cam  = 0;
    
    X = xDim;
    Y = yDim;
        
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

end
