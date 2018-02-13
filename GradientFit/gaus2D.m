function [G] = gaus2D(pos, sig, xid,yid,maxCount)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A = maxCount;
[x,y] = meshgrid(xid,yid);
sigX = sig(1);
sigY = sig(2);
x0 = pos(1);
y0 = pos(2);

xPart = ((x-x0).^2) ./ (2*sigX^2);
yPart = ((y-y0).^2) ./ (2*sigY^2);


G = A.*exp( -(xPart + yPart));

end

