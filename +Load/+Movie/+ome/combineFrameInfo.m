function [frameInfo,totFrame] = combineFrameInfo(frameInfo,fullCombine)
    assert(iscell(frameInfo),'frameInfo is expected to be some cells');
    assert(isstruct(frameInfo{1}),'frameInfo cells are expected to contained struct');
    totFrame = [0 0];
    tmp = [];
    for i = 1 : size(frameInfo,1)
        cFrameInfo = frameInfo{i};
        tmp = [tmp, cFrameInfo];
        totFrame = [totFrame(1) + str2double(cFrameInfo(end-1).T)+1,totFrame(2) + str2double(cFrameInfo(end).T)+1];
    end

    frameInfo = tmp;
    
    if fullCombine
        for i = 1 :size(frameInfo,2)
            if i <=2

            else
                currentT    = str2double(frameInfo(i).T);
                previousT   = str2double(frameInfo(i-1).T);
                previous2T  = str2double(frameInfo(i-2).T);
                currentZ    = str2double(frameInfo(i).Z);
                previousZ   = str2double(frameInfo(i-1).Z);
                currentIFD  = str2double(frameInfo(i).IFD);
                previousIFD = str2double(frameInfo(i-1).IFD);

                %test T
                if currentT < previousT
                    if previous2T == previousT
                        newT = previousT+1;

                    else
                        newT = previousT;
                    end
                    frameInfo(i).T = num2str(newT);
                end

                if currentZ < previousZ

                    newZ = currentZ+1;
                    frameInfo(i).T = num2str(newZ);

                end

                if currentIFD < previousIFD
                    newIFD = previousIFD+1;
                    frameInfo(i).IFD = num2str(newIFD);
                end

            end



        end
    end
    %disp('combining done');
end