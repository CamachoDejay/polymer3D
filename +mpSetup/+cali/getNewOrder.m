function [ neworder, inFocus ] = getNewOrder( inFocus )
%GETNEWORDER get the correct order to place consecutive planes together
%   Detailed explanation goes here

    ordZpos = sort([inFocus.zpos]);
    for i = 1:8
        inFocus([inFocus.zpos] == ordZpos(i)).globalch = i;        
    end
    
    focus = cat(1,inFocus.frame);
    [~,neworder] = sort(focus);
    
end

