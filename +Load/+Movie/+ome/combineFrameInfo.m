function [frameInfo] = combineFrameInfo(frameInfo)
    assert(iscell(frameInfo),'frameInfo is expected to be some cells');
    assert(isstruct(frameInfo{1}),'frameInfo cells are expected to contained struct');
    
    tmp = [];
    for i = 1 : size(frameInfo,1)
        
        tmp = [tmp, frameInfo{i}];
        
    end

    frameInfo = tmp;


end