function [ vals, names ] = rectangulariy( contour )
%Rectangularity represents how rectangular a shape is, i.e. how much it
%fills its minimun biunding rectangle. 
%   Detailed explanation goes here

names = {'Rectangularity'};
vals = cell(length(names),1);

if isempty(contour)
    return
else
    [~,~,Area_shape] = centroid_by_area( contour );
    Area_shape = abs(Area_shape);
    [~, Area_rectangle]= minBoundingBox(contour);

    rec = Area_shape / Area_rectangle;

    vals(1) = {rec};
end

end

