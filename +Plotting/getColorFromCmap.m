%%%%%%%%%%%%%%%%%%%%%%%%%%%% GETCOLORFROMCMAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: get a linearly spaced set of color containing n colors based on
% an existing (char, e.g. jet,hot,...) or inputed colormap (as n x 3
% RGB matrix). This is aimed to be used for plotting multiple graphs on the
% same plot having always differen colors that all belong to a consistent
% colormap.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- cMap - a colormap provided either as a char (matlab existing one,
% see 'doc colormap' for an exhaustive list) or a list of RGB colors in the
% format n x3,
%
% -- n    - the number of linearly spaced color needed 

function colors = getColorFromCmap(cMap,n)

    assert(or(ischar(cMap),ismatrix(cMap)),'colormap needs to be provided as a char(matlab default) or as a matrix of n x 3')
    assert(n==round(n),'n needs to be an integer');

    if ischar(cMap)

        colorArray = colormap(cMap);

    elseif ismatrix(cMap)
        assert(size(cMap,2)==3,'colormap needs to be provided as a char(matlab default) or as a matrix of n x 3')
        colorArray = cMap;

    end

    if n < size(colorArray,1)

        idx = round(linspace(1,size(colorArray,1),n));

        colors = colorArray(idx,:);

    elseif n >= size(colorArray,1)
        
        nIdx = ceil(n/size(colorArray,1));
        nColorArray = zeros(size(colorArray,1)*nIdx,size(colorArray,2));
        for i =1:nIdx
           nColorArray((i-1)*size(colorArray,1)+1:i*size(colorArray)) = ...
               colorArray; 
        end
        
        colors = nColorArray(1:n,:);
        
    end
    
    assert(and(size(colors,1) == n, size(colors,2) ==3),...
        'Unexpected output dimension, something went wrong');
       
end