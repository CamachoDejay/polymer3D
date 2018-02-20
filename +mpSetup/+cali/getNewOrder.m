function [ neworder, in_focus ] = getNewOrder( in_focus )
%GETNEWORDER get the correct order to place consecutive planes together
%   Detailed explanation goes here

    or_zpos = sort([in_focus.zpos]);
    for i = 1:8
        in_focus([in_focus.zpos] == or_zpos(i)).globalch = i;        
    end
    
    neworder = zeros(8,1);
    for i = 1:8
        neworder(i) = find([in_focus.globalch]==i);
    end

end

