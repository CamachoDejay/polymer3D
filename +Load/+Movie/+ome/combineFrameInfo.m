function [frameInfo,totFrame] = combineFrameInfo(frameInfo)
    assert(iscell(frameInfo),'frameInfo is expected to be some cells');
    assert(isstruct(frameInfo{1}),'frameInfo cells are expected to contained struct');
    totFrame = [0 0];
    tmp = [];
    %pilling up data from different files on top of each other, sorting and counting
    %total number of frame
    for i = 1 : size(frameInfo,1)
        cFrameInfo = frameInfo{i};
        fCam1 = {cFrameInfo(:,contains({cFrameInfo.C},'0')).T};
        fCam2 = {cFrameInfo(:,contains({cFrameInfo.C},'1')).T};
        maxFr = length(cFrameInfo)/2;
        %fix if there is only one frame in the frame
        if and(length(unique([cFrameInfo.T])) == 1,length(unique([cFrameInfo.Z])) == 1)
            
            if i>1
                %Change data depending on previous
                nFrames = str2double(prevFInfo(1).T)+1;
                cFrameInfo(1).T = num2str(nFrames);
                cFrameInfo(2).T = num2str(nFrames);
                
                nFrames = str2double(prevFInfo(1).Z)+1;
                cFrameInfo(1).Z = num2str(nFrames);
                cFrameInfo(2).Z = num2str(nFrames);
                
            end
                  
            prevFInfo = cFrameInfo;
            %need to modify cFrameInfo so correct data is saved.
        else
            data = cellfun(@str2double,{cFrameInfo.T});
            
            [dataSorted,sortIdx] = sort(data);
            
            cFrameInfo = cFrameInfo(sortIdx);
            if i>1
                if cFrameInfo(1).T < tmp(end).T
                    dataSorted = dataSorted + (i-1)*(maxFr);
                    dataSorted = num2cell(dataSorted);
                    dataSorted = cellfun(@num2str,dataSorted,'UniformOutput',false);
                    for j = 1: length(dataSorted)
                        cFrameInfo(j).T = dataSorted{j};
                    end
                    
                    
                end
            end
            
            
        end
        tmp = [tmp, cFrameInfo];
        
        
        totFrame = [totFrame(1) + str2double(fCam1{end})+1,totFrame(2) + str2double(fCam2{end})+1];
    end
    
    frameInfo = tmp;
    
    assert(totFrame(1) == totFrame(2),'Number of frame of both camera are not consistent');
    %disp('combining done');
end