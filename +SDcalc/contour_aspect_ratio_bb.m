function [ vals, names ] = contour_aspect_ratio_bb( contour )
% Calculation of the aspect ratio using the minimun bounding box method.
%   Input:  2xn matrix containing the [x,y] coordinates of n points

names = {'Aspect_Ratio'; 'Max_Feret'; 'Min_Feret'};
vals = cell(length(names),1);

if isempty(contour)
    return
else


    assert(size(contour,1)==2,['Input must be a 2xn matrix containing'...
                               ' the [x,y] coordinates of n points'])
    assert(size(contour,2)>3,['Input must be a 2xn matrix containing'...
                               'the [x,y] coordinates of n points and a'...
                        ' bounding box cant be defined for less than 3 points'])
    [bBox, ~]= minBoundingBox(contour);
    E = diff(bBox,1,2);
    d = sqrt(  sum(E.^2,1) );
    L = max(d);
    W = min(d);
    aspect_ratio = L/W;
    
    vals(1) = {aspect_ratio};
    vals(2) = {L};
    vals(3) = {W};
end

end

