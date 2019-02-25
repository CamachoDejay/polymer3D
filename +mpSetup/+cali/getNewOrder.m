function [ neworder, inFocus ] = getNewOrder( inFocus )
%GETNEWORDER get the correct order to place consecutive planes together
%   Detailed explanation goes here

    ordZpos = sort([inFocus.zpos],'descend');%Descending order because the 
    %first plane to get into focus is the plane that is the highest in Z.
    for i = 1:length(inFocus)
        inFocus([inFocus.zpos] == ordZpos(i)).globalch = i;        
    end
    
    focus = cat(1,inFocus.frame);
    [~,neworder] = sort(focus,'descend');
    
end

