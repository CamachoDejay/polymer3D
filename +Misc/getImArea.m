function [out] = getImArea(IM)
%GETPREPOST Summary of this function goes here
%   Detailed explanation goes here
% [info] = loadMovie.tif.getinfo(fPath);

figure(1)
clf
imagesc(IM)
title('Full image')
axis image
drawnow;
waitfor(helpdlg('Choose area to crop'));
rect = getrect;
rect = round(rect);
xi = rect(2);
xf = xi + rect(4);
yi = rect(1);
yf = yi + rect(3);
out = IM(xi:xf,yi:yf);

end

