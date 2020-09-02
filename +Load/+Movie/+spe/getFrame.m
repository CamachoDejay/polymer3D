%% Importing data from .SPE file into matlab array 'Im' in parts
%
% Allows for frame specific import of SPE files when the entire file may
% either be too large to handle or if only a specific part of the image
% stack is important. Can be used in a for loop to divide a large stack of
% images into multiple arrays of either the same or different sizes.
%
% Ex:
%
% Im=ImportSPEframes(Ni,Nf);
%
% Ni is the starting frame and Nf is the ending frame of image stack Im. An
% error is returned if the ending frame request is larger than the actual
% image stack. Setting Ni=1 and Nf=0 returns the entire image stack.
%


function Im = getFrame(filePath,frames)
    assert(ischar(filePath),'Error, wrong path Names');
    assert(min(size(frames)) ==1, 'Frames to import should be a vector');
    
    fid = fopen(filePath,'r'); 

    header = fread(fid,2050,'uint16=>uint16'); % 2050 uint16 = 4100 bytes = 32800 bits

    Xdim = double(header(22));
    Ydim = double(header(329));
    Zdim = double(header(724));
    DataType = header(55);
    
    switch nargin
        case 1
             frames = 1:Zdim;
             nF = length(frames);
             nI = 1;
             disp('Loading the entire stack of images, let us hope we have enough memory');
             
        case 2
            nF = frames(end);
            nI = frames(1);
    end
    
    assert(max(frames)<=Zdim,'Requested number of frame exceeds the file');
    
    switch DataType
        case 0	% FLOATING POINT (4 bytes / 32 bits)
            ImMat = fread(fid,inf,'float32=>float32');
        case 1	% LONG INTEGER (4 bytes / 32 bits)
            ImMat = fread(fid,inf,'int32=>int32');
        case 2	% INTEGER (2 bytes / 16 bits)
            ImMat = fread(fid,inf,'int16=>int16');
        case 3	% UNSIGNED INTEGER (2 bytes / 16 bits)
            if nI~=1
                SKIP=(nI-1)*Xdim*Ydim;
                fseek(fid,2*(SKIP+2050),'bof');
                
                Range=(nF-nI+1)*Xdim*Ydim;
                
                ImMat = fread(fid,Range,'uint16=>uint16');
                Zdim=size(ImMat,1)/(Xdim*Ydim);
            else
                Range=(nF-nI+1)*Xdim*Ydim;
                ImMat = fread(fid,Range,'uint16=>uint16');
                Zdim=size(ImMat,1)/(Xdim*Ydim);
            end
            
    end

    fclose(fid);

    a = reshape(ImMat,Xdim,Ydim,Zdim);
    % clear ImMat; %clear some memory

    %permute the X and Y dimensions so that an image looks like in Winview
    a = permute(a,[2,1,3]);
    Im = double(a);
end