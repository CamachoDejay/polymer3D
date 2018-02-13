% extract the candidate ROI
% -------------------------------------------------------------------------
% Inputs
%   imRAW      : Raw image
%   imDN       : Denoised image
%   RegR       : Region radius of the ROI, so the size of ROI is (2*WinR)+1
%   threshold  : Theshold to pick up the candidate emitter
% -------------------------------------------------------------------------
% Outputs
%   ROIs       : The cadidate emitters stack
%   ROI_coor   : The coordinates of the cadidate ROIs (Pixels) 
% -------------------------------------------------------------------------

function [ROIs, ROI_coor] = extractROI(imRAW,imDN,RegR,threshold)

[imW,imH] = size(imDN);
ROIs = zeros(2*RegR+1,2*RegR+1,1);
ROI_coor = zeros(2,1);
ROInum = 0;
for m=1+RegR:imW-RegR-1
    for n=1+RegR:imH-RegR-1
        if (imDN(m,n)>threshold) && (imDN(m,n)==max(max(imDN(m-2:m+2,n-2:n+2))))
            imDN(m,n) = imDN(m,n) + 1e-6;
            ROInum = ROInum + 1;
            ROIs(:,:,ROInum) = imRAW(m-RegR:m+RegR,n-RegR:n+RegR);
            ROI_coor(:,ROInum) = [m,n];
        end
    end
end
       
end