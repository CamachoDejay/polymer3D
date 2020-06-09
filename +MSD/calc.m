function [MSD,D] = calc(coord)
% the function gives the distance ('D') between each point and
% mean-squre-dispracement ('MSD')
%
%input parameters; cod is 2D or 3D aray including coordinate of trajectry.

dim = size(coord,2);

switch dim
    case 1
        
        coord(:,2:3) = 0;
        
    case 2
        
        coord(:,3) = 0;
        
    case 3
        
        
    otherwise
        
        error('unexpected dimension for the vector')
        
end

%%%Fist calculate distance between each frame

DX = diff(coord(:,1).');
DY = diff(coord(:,2).');
DZ = diff(coord(:,3).');
D = sqrt(DX.^2 + DY.^2 + DZ.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MSD = zeros(size(coord,1)-1,1);
%Calculate mean-squere-displacement
for i = 1:size(coord,1)-1
    
    stp = i;
    cnt =  1;
    D1  = [];
    while cnt<=stp && cnt+stp<=size(coord,1)
        
        idx = cnt:stp:size(coord,1);
        DX  = diff(coord(idx,1).');
        DY  = diff(coord(idx,2).');
        DZ  = diff(coord(idx,3).');
        D1  = [D1 sqrt(DX.^2 + DY.^2 + DZ.^2)];
        cnt = cnt+1;
        
        if ~isempty(D1)
            
            D2=D1(~isnan(D1));
            
            if ~isempty(D2)
                
                MSD(i) = mean(D2.^2);
                
            else
                
                MSD(i) = NaN;
                
            end
        end
    end %while
end

D   = D(:);
MSD = MSD(:);

end