function [BW] = globThresh(imStack,connectivity,diskDim)
        
    BW = imbinarize(imStack,'global');
    BW = ~BW;
    BW = bwareaopen(BW,connectivity);
    SE = strel('disk',diskDim);
    BW = imopen(BW,SE);
    
end