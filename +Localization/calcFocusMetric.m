function [corrEllip, focusMetric] = calcFocusMetric(ellipticity,LRT,weight)
assert(isvector(ellipticity),'Ellipticity should be given as a vector');
switch nargin
    
    case 1 
        
        LikelyH = false;
        
    case 2
        
        LikelyH = true;
        weight = 1;
        
    case 3
        
        assert(and(isnumeric(weight),max(size(weight))==1), 'Weight should contain only one number')
        LikelyH = true;
        
end
corrEllip = ellipticity;
corrEllip(corrEllip>1)= 1./ellipticity(ellipticity>1);

if LikelyH
    assert(isvector(LRT),'Likelyhood ratio test should be given as a vector');
    assert(length(ellipticity) == length(LRT),'Input vector should have the same size');
    
    focusMetric = abs(LRT) .* corrEllip.^weight;
    
else
    
    focusMetric = [];
    
end

end