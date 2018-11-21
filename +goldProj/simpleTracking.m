function [newOrder] = simpleTracking(newPos,prevPos)

nPart = length(newPos)/2;
newOrder = zeros(nPart,1);
minVal = zeros(nPart,1);
euclDist = zeros(nPart);

for i = 1 : nPart
    newY = newPos(2*i);
    newX = newPos(2*i-1);
    prevY = prevPos(2:2:end);
    prevX = prevPos(1:2:end);
    
    euclDist(i,:) = sqrt((newX-prevX).^2 +(newY-prevY).^2);
    
%     [ minVal(i),newOrder(i)] = min(euclDist(i,:));
end

[newOrder,~] = goldProj.munkres(euclDist');
% u=unique(newOrder);         % the unique values
% [n,bin]=histc(newOrder,u);  % count how many of each and where
% ix1 = find(n>1,1); 
% 
% if ~isempty(ix1)
% disp('test')
% newOrder = zeros(nPart,1);
% 
% 
% 
% 
% end
% 

end