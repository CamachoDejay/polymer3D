function [gSegStack,aSegStack] = segmentStack(imStack,varargin)
% segment an 2D or 3D image.
%
% function [gSegStack,aSegStack] = segmentStack(imStack, ...)
%
% Purpose 
% Segment images for further analysis (e.g. polymer network)
%
%
%    Inputs (required)
% imStack - 2D or 3D matrix where the two first dimension represent an
% image(x,y). 3rd dimension is variation of the image in space (z)
%
%   Inputs (optional, param/value pairs)
% 'frames' - [imSize(3) by default]define frame range to be analyze
%
% 'connectivity' -[8 or 216 by default for 2D and 3D segmentation respectively
% define the pixel connectivity needed to belong to the same cluster?
% 
% 'method' - ['both' by default]: 'global', 'adaptive' or 'both'
% segmentation method
%
% 'threshold' - [0.5 by default] defines the sensitivity for adaptive
% threshold
%          
% 'diskDim' - [4  by default] defines the diameter of the disk for cleaning
% up the binary image

% 'neigh'  - [301 301 151 by default] the number of neighbor to be consider
% for adaptive threshold

%   Outputs
%'gSegStack' - binarized image with global method
%'aSegStack' - binarized image with adaptive method


%parse user input
narginchk(1,inf)

%Calculate size and dimensions of the image stack
imSize = size(imStack);
dim = length(imSize);
imSize(3) = size(imStack,3);

assert(ismember(dim,[2 3]),'imStack is expected to be a 2D or 3D matrix');

params = inputParser;
params.CaseSensitive = false;
params.addParameter('frames', 1:imSize(3), @(x) isnumeric(x) && isvector(x) && all(x>0));
params.addParameter('connectivity', 8*double(dim==2)+216*double(dim==3),...
    @(x) isnumeric(x) || x> 0);
params.addParameter('method', 'both',...
    @(x) ischar(x) && (strcmp(x,'both') || strcmp(x,'adaptive') || strcmp(x,'global')));
params.addParameter('threshold',0.5,@(x) isnumeric(x) && x<=1 && x>0);
params.addParameter('diskDim', 4, @(x) isnumeric(x) && x>1 );
params.addParameter('neigh',[301 301 151],@(x) all(isnumeric(x)) && length(size(x))==dim && all(x>1));

params.parse(varargin{:});

%Extract values from the inputParser
frames =  params.Results.frames;
connectivity =  params.Results.connectivity;
method = params.Results.method;
threshold = params.Results.threshold;
diskDim = params.Results.diskDim;
neigh = params.Results.neigh;

imStack = imStack(:,:,frames);

%act depending on the method chosen by the user
switch method
    case 'global'
        
        [BWglobal] = globThresh(imStack,connectivity,diskDim);
        gSegStack = BWglobal;
        aSegStack = [];
    case 'adapt'
        
        [BWadapt] = adaptiveThresh(imStack,connectivity,threshold,diskDim,neigh);
        gSegStack = [];
        aSegStack = BWadapt;
        
    case 'both'
        
        [BWglobal] = globThresh(imStack,connectivity,diskDim);
        [BWadapt]  = adaptiveThresh(imStack,connectivity,threshold,diskDim,neigh);
        aSegStack  = BWadapt;
        gSegStack  = BWglobal;
    otherwise
        
        error('unknown segmentation method requested, only know "global", "adaptive", "both"');

end

end


function [BW] = globThresh(imStack,connectivity,diskDim)
        
    BW = imbinarize(imStack,'global');
    BW = ~BW;
    BW = bwareaopen(BW,connectivity);
    SE = strel('disk',diskDim);
    BW = imopen(BW,SE);
    
end

function [BW] = adaptiveThresh(imStack,connectivity,threshold,diskDim,neigh)

        if size(imStack,3) ==1
            th = adaptthresh(imStack,threshold,'Fore','bright');
        else
            th = adaptthresh(imStack,threshold,'neigh',neigh,'Fore','bright');
        end
        
        BW = imbinarize(imStack,th);
        BW = ~BW;
        BW = bwareaopen(BW,connectivity);
        SE = strel('disk',diskDim);
        BW = imopen(BW,SE);
        
end