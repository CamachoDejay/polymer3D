function [ data12 ] = arrange_correctInt( data12, cal )
%ARRANGE_CORRECTINT rearranges data and corrects for intensity diffecences
%between channels
%   Detailed explanation goes here

    C     = cal.Icorrf;
    newor = [cal.neworder];
    % correct int
    for i = 1:8
        data12(:,:,i,:) = data12(:,:,i,:)./C(i);
    end
    
    % reorder planes
    data12 = data12(:,:,newor,:);
    
end

