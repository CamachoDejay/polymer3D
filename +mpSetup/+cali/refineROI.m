function [ ROI ] = refineROI( ROI, im_shifts )
%REFINEROI improves the ROI to have the smallest shifts possibles
%   Detailed explanation goes here
    
    row_w = ROI(1,4);
    col_w = ROI(1,3);
    %get min and max imShift for row
    mars = max(im_shifts(:,1));
    mirs = min(im_shifts(:,1));
    %get in and max imShift for Col
    macs = max(im_shifts(:,2));
    mics = min(im_shifts(:,2));
    %calculate the width 
    row_w = row_w - mars + mirs;
    col_w = col_w - macs + mics;
    
    row_i = ROI(:,2)+(mars-(im_shifts(:,1)));
    col_i = ROI(:,1)+(macs-(im_shifts(:,2)));
    
    ROI(:,4) = row_w;
    ROI(:,2) = row_i;
    ROI(:,3) = col_w;
    ROI(:,1) = col_i;

end

