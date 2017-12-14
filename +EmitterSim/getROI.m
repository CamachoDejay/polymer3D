function [roi_lims] = getROI(x0_pix, y0_pix, delta, x_size, y_size)
%GETROI Summary of this function goes here
%   Detailed explanation goes here

x1 = round(x0_pix - delta);
x2 = round(x0_pix + delta);
y1 = round(y0_pix - delta);
y2 = round(y0_pix + delta);

if x1 < 1;
    x1 = 1;
end

if x2 > x_size
    x2 = x_size;
end

if y1 < 1
    y1 = 1;
end

if y2 > y_size
    y2 = y_size;
end

roi_lims = [x1 x2 y1 y2];

end

