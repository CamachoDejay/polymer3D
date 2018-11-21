function [pos] = simpleRegMaxDetection (data,nMax)
    
    x = 1:size(data,2);
    y = 1:size(data,1);
    [domX,domY] = meshgrid(x,y);

    %rough estimate of bkg:
    bkg = median(median(data));
    data = data-bkg; %bkg subtraction
    data(data<0.1*max(max(data))) = 0;
    out = imregionalmax(data);

    x0 = domX(out);
    x0 = x0(x0~=0);
    y0 = domY(out);
    y0 = y0(y0~=0);
    bright = data(out);
    bright = bright(bright~=0);
    %find the brightess maxima according to the number of particles specified
    [~,idx] = maxk(bright,nMax);
    x0 = x0(idx);
    y0 = y0(idx);
    
    pos = [y0 x0];

end