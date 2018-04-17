function [ROIs] = getROIs(locList,ROIrad,imSize)
%GETROIS get ROIs from list of localizations
%   Detailed explanation goes here

% build ROIs
ROIcent = round(locList);
nLoc   = size(ROIcent,1);
% imSize = [size(data,1),size(data,2)];

% test for molecules to close to the im edges
for i = 1:2
    % test for loc that are too close to the beginning of the image
    tVal = (ROIcent(:,i) - ROIrad);
    test = tVal < 1;
    ROIcent(test,i) = ROIrad +1;
    
    % test for loc that are too cose to the end of the im
    tVal = (ROIcent(:,i) + ROIrad);
    test = tVal > imSize(i);
    ROIcent(test,i) = imSize(i) - ROIrad;
    
end

% storing of ROIs
ROIs = zeros(nLoc,4);
ROIs(:,1) = round(ROIcent(:,1))-ROIrad;
ROIs(:,2) = round(ROIcent(:,1))+ROIrad;
ROIs(:,3) = round(ROIcent(:,2))-ROIrad;
ROIs(:,4) = round(ROIcent(:,2))+ROIrad;

imLims = repmat([1 imSize(1) 1 imSize(2)],nLoc,1);
test = [ROIs(:,1)>=imLims(:,1), ROIs(:,2)<=imLims(:,2),...
        ROIs(:,3)>=imLims(:,3), ROIs(:,4)<=imLims(:,4)];
    
assert(all(test(:)), 'Problems, unexpected issues during ROI definition')

end

