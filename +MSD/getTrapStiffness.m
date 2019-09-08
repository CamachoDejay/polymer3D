function [trapStiffness] = getTrapStiffness(coord,T)
    Kb = 1.381e-23;
    assert(size(coord,1)>size(coord,2), 'Coordinate needs to be provided as column vector');
    coord = coord - mean(coord,1);
    
    displacement =  mean(diff(coord*10^(-9),1).^2,1);
    
    % to convert to pN/um we need to multiply by 10e-12 and divide by 10e-6
    % 
    trapStiffness = Kb*T./displacement/10^-6;

end