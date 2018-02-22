function [ c1, c2 ] = consolidatePos( mainSet, testSet, rTh )
% consolidates positions that are at a distance smaller than rTh
%   Detailed explanation goes here

[IDX, D] = knnsearch(mainSet,testSet);
common = (D<rTh);


c1   = mainSet(IDX(common),:);
c2   = testSet(common,:);

end

