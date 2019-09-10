% This function aim to calculate the trap stiffness based on temperature
% and variance in displacement. The equation can be found here:
% http://genomics.princeton.edu/shaevitzlab/OT_Practicle_Guide.pdf
% 
% INPUT
%   -  coord: set of coordinate (1,2 or 3 column for 1,2 or 3 dimension)
%             positions need to be in meter beforehand
%   -  T: Temperature in Kelvin
%         
% OUTPUT
%   -  trapStiffness: Stiffness of the trap (spring constant) in pN/um;

function [trapStiffness] = getTrapStiffness(coord,T)
    assert(size(coord,1)>size(coord,2), 'Coordinate needs to be provided as column vector');    
    
    Kb = 1.381e-23;
    coord = coord - mean(coord,1);
    var          = (sum(coord.^2,1))./size(coord(1:end-1,1),1);
    
    % to convert to pN/um we need to multiply by 10e-12 and divide by 10e-6
    trapStiffness = Kb*T./var/10^-6;

end