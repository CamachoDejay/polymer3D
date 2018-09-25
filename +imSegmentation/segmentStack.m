function [gSegStack,aSegStack] = segmentStack(imStack,frames,connectivity,method,threshold,diskDim)
%Calculate size and dimensions of the image stack
imSize = size(imStack);
dim = length(imSize);
assert(ismember(dim,[2 3]),'imStack is expected to be a 2D or 3D matrix');

%Check number of input argument and provide default input if not provided
switch nargin
    case 1
        frames = 1:size(imStack,3);
        threshold = 0.5;
        method = 'both';
        if dim == 3
            connectivity = 216;
        else
            connectivity = 8;
        end
        diskDim = 4;
    case 2
        
        threshold = 0.5;
        method = 'both';
        if dim == 3
            connectivity = 216;
        else
            connectivity = 8;
        end
        diskDim = 4;
    case 3
        
        threshold = 0.5;
        method = 'both';
        diskDim = 4;
        
    case 4
        
        threshold = 0.5;
        diskDim = 4;
        
    case 5
        
        diskDim = 4;
    case 6
        
    otherwise
        
        error('Too many input arguments')
        
end
%Check User input            
assert(max(frames) <= size(imStack,3),'requested number of frame exceed max frame');
assert(isnumeric(connectivity),'connectivity is expected to be a number');
assert(ischar(method),'method is expected to be a chain of char, valid input are: global, adaptive, both');
assert (isnumeric(threshold),'threshold is expected to be a number');
assert(and(threshold >=0, threshold <= 1),'threshold is expected to be comprised between 0 and 1');

%act depending on the method chosen by the user
switch method
    case 'global'
        
        [BWglobal] = imSegmentation.globThresh(imStack,connectivity,diskDim);
        gSegStack = BWglobal;
        aSegStack = [];
    case 'adapt'
        
        [BWadapt] = imSegmentation.adaptiveThresh(imStack,connectivity,threshold,diskDim);
        gSegStack = [];
        aSegStack = BWadapt;
        
    case 'both'
        
        [BWglobal] = imSegmentation.globThresh(imStack,connectivity,diskDim);
        [BWadapt]  = imSegmentation.adaptiveThresh(imStack,connectivity,threshold,diskDim);
        aSegStack  = BWadapt;
        gSegStack  = BWglobal;
    otherwise
        
        error('unknown segmentation method requested, only know "global", "adaptive", "both"');

end

end
