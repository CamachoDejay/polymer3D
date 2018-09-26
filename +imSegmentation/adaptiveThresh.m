function [BW] = adaptiveThresh(imStack,connectivity,threshold,diskDim)
        
        dim = length(size(imStack));
        
        switch dim
            case 2
                th = adaptthresh(imStack,threshold,'Fore','bright');
            case 3
                th = adaptthresh(imStack,threshold,'neigh',[301 301 151],'Fore','bright');
            otherwise
                error('dimension of the image to segment does not make sense');
        end
        
        BW = imbinarize(imStack,th);
        BW = ~BW;
        BW = bwareaopen(BW,connectivity);
        SE = strel('disk',diskDim);
        BW = imopen(BW,SE);
        
end