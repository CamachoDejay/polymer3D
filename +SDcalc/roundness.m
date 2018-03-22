function [ vals, names ] = roundness( contour )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

names = {'Roundness_Perim';'Roundness_Area'};
vals = cell(length(names),1);

if isempty(contour)
    return
else
    [Perim]    = perimeter_li( contour );
    Perim = Perim{1};
    [~,~,Area] = centroid_by_area( contour );
    Area       = abs(Area);

    metric_perim = Perim / (2*(pi*Area)^0.5);
    metric_area  = 4 * pi * Area / ( Perim^2 );
    
    vals(1) = {metric_perim};
    vals(2) = {metric_area};

% keep in mind that:
% metric_perim = 1 / (metric_area)^2
% thus this two metrics are not independent from each other.
end

end

