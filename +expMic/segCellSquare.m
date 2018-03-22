function [C, R] = segCellSquare(IM,tHold,titleStr)
%SEGCELLSQUARE segment the cell and the square
%   Detailed explanation goes here

% apply treshold
BWpre = IM>tHold;
% clear border
BWpre = imclearborder(BWpre);
% remove small things
BWpre = bwareaopen(BWpre,5);
% smooth
se = strel('disk',3);
BWpre = imclose(BWpre,se);
BWpre = imopen(BWpre,se);
% remove small things
BWpre = bwareaopen(BWpre,100);
BWpre = ~BWpre;
BWpre = bwareaopen(BWpre,100);
BWpre = ~BWpre;

% get boundaries
[B,~] = bwboundaries(BWpre);
% size of the boundries
Blength = cellfun(@length,B);
% cell and square must be the largest of them
[~, idx] = maxk(Blength,2);
B = B(idx);
% Blength = Blength(idx);
assert(length(B)==2,'We seem to find more than 2 contours, we dont handle that at the moment but it can be done')

%resample and smooth
np = 500;
B{1} = SDcalc.smoothBoundary( B{1} );
B{1} = SDcalc.resampleBoundary( B{1}, np );
B{1} = SDcalc.smoothBoundary( B{1} );

B{2} = SDcalc.smoothBoundary( B{2} );
B{2} = SDcalc.resampleBoundary( B{2}, np );
B{2} = SDcalc.smoothBoundary( B{2} );

% looking at the first area and getting shape descriptors
contour = B{1}';
rec     = SDcalc.rectangulariy(contour);
[x,y,Ar]= SDcalc.centroid_by_area( contour);
[bb, ~] = SDcalc.minBoundingBox(contour);
E       = diff(bb,1,2);
d       = sqrt(  sum(E.^2,1) );
d       = d(1:2);
[L,i]   = max(d);
line    = bb(:,i:i+1);
dl      = diff(line,1,2);
ang     = atan(dl(2)/dl(1))*180/pi;
ang     = mod(ang,(180));

ROI1.rec = rec{1};
ROI1.contour = contour;
ROI1.centroid = [x,y];
ROI1.area = abs(Ar);
ROI1.boundingBox = bb;
ROI1.length = L;
ROI1.ang = ang;
ROI1.contNorm = contour-repmat([x;y],1,np+1);

% looking at the second area and getting shape descriptors
contour = B{2}';
rec     = SDcalc.rectangulariy(contour);
[x,y,Ar]= SDcalc.centroid_by_area( contour);
[bb, ~] = SDcalc.minBoundingBox(contour);
E       = diff(bb,1,2);
d       = sqrt(  sum(E.^2,1) );
d       = d(1:2);
[L,i]   = max(d);
line    = bb(:,i:i+1);
dl      = diff(line,1,2);
ang     = atan(dl(2)/dl(1))*180/pi;
ang     = mod(ang,(180));

ROI2.rec = rec{1};
ROI2.contour = contour;
ROI2.centroid = [x,y];
ROI2.area = abs(Ar);
ROI2.boundingBox = bb;
ROI2.length = L;
ROI2.ang = ang;
ROI2.contNorm = contour-repmat([x;y],1,np+1);

if ROI2.area<ROI1.area
    % then the second area is the square
    R = ROI2;
    C = ROI1;
else
    % then the first area is the square
    R = ROI1;
    C = ROI2;
end

if C.rec > R.rec
    warning('cell looks more squared than square')
end

% image 
figure(1)
clf
imagesc(IM)
axis image
title(titleStr)

hold on
plot(R.contour(2,:), R.contour(1,:), ':r', 'LineWidth', 2)
plot(C.contour(2,:), C.contour(1,:), ':g', 'LineWidth', 2)
hold off
end

