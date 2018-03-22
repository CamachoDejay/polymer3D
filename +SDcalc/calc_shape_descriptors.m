function [ vals, names ] = calc_shape_descriptors (contour)
%This function calculates the descriptors for a given contour
%   Detailed explanation goes here

    vals = [];
    names = [];
    [ pval, pname ] = perimeter_li( contour );
    vals = [vals; pval];
    names = [names; pname];

    [ rovals, ronames ] = roundness( contour );
    vals = [vals; rovals];
    names = [names; ronames];
    
    [ apvals, arnames ] = contour_aspect_ratio_bb( contour );
    vals = [vals; apvals];
    names = [names; arnames];

    [ revals, renames ] = rectangulariy( contour );
    vals = [vals; revals];
    names = [names; renames];
    
    [ cvals, cnames ] = convexity( contour );
    vals = [vals; cvals];
    names = [names; cnames];
    
    [ svals, snames ] = solidity( contour );
    vals = [vals; svals];
    names = [names; snames];
        
    [ ENvals, ENnames ] = EnergyDescriptors( contour );
    vals = [vals; ENvals];
    names = [names; ENnames];
    
    
    [c_dist] = centroid_distance( contour );
    r_in      = c_dist.fit.r;
    t_in      = c_dist.fit.t;
    n_comp    = 20;
    [ FD, FDnames ] = FD_centroid_distance( r_in, t_in, n_comp, false );
    if isempty(FD.components)
        FDvals = cell(length(FDnames),1);
    else
        FDvals = num2cell(FD.components);
    end
    vals = [vals; FDvals];
    names = [names; FDnames];
 
    
end

