function [pos] = simpleRegMaxDetection (data,nMax,minDist)
    
    x = 1:size(data,2);
    y = 1:size(data,1);
    [domX,domY] = meshgrid(x,y);

    %rough estimate of bkg:
    bkg = median(median(data));
    data = data-bkg; %bkg subtraction
    data(data<0.1*max(max(data))) = 0;
    %gaussian filtering
   % data = imgaussfilt(data,3);
    
    out = imregionalmax(data);

    x0 = domX(out);
    x0 = x0(x0~=0);
    y0 = domY(out);
    y0 = y0(y0~=0);
    bright = data(out);
    bright = bright(bright~=0);
    
    %Remove maxima that are too close to each others
    for i = 1:length(x0)
        
        if i > length(x0)
            break;
        end
        
        euclDist = sqrt((x0(i) -x0).^2 + (y0(i) -y0).^2);
        idx = find(euclDist<minDist);
        
        %check which is brightest
        [~,idx2] = max(bright(idx));
        
        %delete them
        x0(bright(idx) < bright(idx(idx2))) = [];
        y0(bright(idx) < bright(idx(idx2))) = [];
        bright(bright(idx) < bright(idx(idx2))) = [];

    end
    
    %find the brightess maxima according to the number of particles specified
    [~,idx] = maxk(bright,nMax);
    x0 = x0(idx);
    y0 = y0(idx);
    
    pos = [y0 x0];

end