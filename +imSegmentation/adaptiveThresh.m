function [BW] = adaptiveThresh(imStack,connectivity,threshold,diskDim)

        th = adaptthresh(imStack,threshold,'neigh',[301 301 151],'Fore','bright');
        BW = imbinarize(imStack,th);
        BW = ~BW;
        BW = bwareaopen(BW,connectivity);
        SE = strel('disk',diskDim);
        BW = imopen(BW,SE);
        
end