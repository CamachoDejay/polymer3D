function [newOrder] = simpleTracking(newPos,prevPos)

nPart = length(newPos)/2;
sqDist = zeros(nPart);

for i = 1 : nPart
    newY = newPos(2*i);
    newX = newPos(2*i-1);
    prevY = prevPos(2:2:end);
    prevX = prevPos(1:2:end);
    sqDist(i,:) = (newX-prevX).^2 +(newY-prevY).^2;
  
end

[newOrder,~] = Minimization.munkres(sqDist');


end