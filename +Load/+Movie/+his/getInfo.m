function [frameInfo,movInfo] = getInfo(path2File)
    
    %Create a bioformat reader
    imReader = bfGetReader(path2File);
    %Get the metaData
    metaData = imReader.getSeriesMetadata();
    
    %get the size of the image
    sizeData = metaData.get('vrcAcqArea');
    
    idx = strfind(sizeData,',');
    %extract dimension
    Z = str2double(sizeData(1:idx(1)-1));
    T = str2double(sizeData(idx(1):idx(2)-1));
    X = str2double(sizeData(idx(2)+1:idx(3)-1));
    Y = str2double(sizeData(idx(3)+1:end));
    
    %get the number of frame
    maxFrame = imReader.getMetadataStore().getImageCount();
    %get the exposure time
    expT     = str2double(metaData.get('vExpTim1'))*1000; %in ms
    %store a few thing for the rest of the code
    isMultiImage = false;
    
    isZStack = false;
    
    Cam  = 0;
    
    %Store info for output
    frameInfo = [];
    
    movInfo.Width  = X;
    movInfo.Length = Y;
    movInfo.Path   = path2File;
    movInfo.isMultiImage = isMultiImage;
    movInfo.isZStack = isZStack;
    movInfo.Cam = Cam;
    movInfo.expT = expT;
    movInfo.maxFrame = maxFrame;
    
    
    
    
    
    


end