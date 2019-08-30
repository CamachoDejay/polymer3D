function [newOrder] = simpleTracking(newPos,prevPos,dim)

switch dim
    
    case '2D'
        nPart = length(newPos)/2;
        sqDist = zeros(nPart);

        for i = 1 : nPart
            newY = newPos(2*i);
            newX = newPos(2*i-1);
            prevY = prevPos(2:2:end);
            prevX = prevPos(1:2:end);
            sqDist(i,:) = (newX-prevX).^2 +(newY-prevY).^2;

        end
    case '3D'
        nPart = length(newPos)/3;
        sqDist = zeros(nPart);

        for i = 1 : nPart
            newZ = newPos(3*i);
            newY = newPos(3*i-1);
            newX = newPos(3*i-2);
            prevZ = prevPos(3:3:end);
            prevY = prevPos(2:3:end);
            prevX = prevPos(1:3:end);
            sqDist(i,:) = (newX-prevX).^2 +(newY-prevY).^2 +(newZ-prevZ).^2;
        end
        
end

[newOrder,~] = Minimization.munkres(sqDist');


end