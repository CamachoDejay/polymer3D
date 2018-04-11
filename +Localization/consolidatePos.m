function [ mainSet ] = consolidatePos( mainSet, testSet, rTh )
% consolidates positions that are at a distance smaller than rTh
%   Detailed explanation goes here

[IDX, D] = knnsearch(mainSet,testSet);
common = (D<rTh);


c1   = mainSet(IDX(common),:);
c2   = testSet(common,:);
ccom = (c1+c2)./2;
cnew = testSet(~common,:);

mainSet(IDX(common),:) = ccom;
mainSet = cat(1,mainSet,cnew);


end

