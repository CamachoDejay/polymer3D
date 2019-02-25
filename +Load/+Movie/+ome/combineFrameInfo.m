function [frameInfo,totFrame] = combineFrameInfo(frameInfo,fullCombine)
    assert(iscell(frameInfo),'frameInfo is expected to be some cells');
    assert(isstruct(frameInfo{1}),'frameInfo cells are expected to contained struct');
    totFrame = [0 0];
    tmp = [];
    for i = 1 : size(frameInfo,1)
        cFrameInfo = frameInfo{i};
        tmp = [tmp, cFrameInfo];
        
        fCam1 = {cFrameInfo(:,contains({cFrameInfo.C},'0')).T};
        fCam2 = {cFrameInfo(:,contains({cFrameInfo.C},'1')).T};
        totFrame = [totFrame(1) + str2double(fCam1{end})+1,totFrame(2) + str2double(fCam2{end})+1];
    end
    
    frameInfo = tmp;
    
    if fullCombine
        
        maxFr = max(cellfun(@str2double,fCam1));
        bothCam = (maxFr+1)*2;
        nIt = size(frameInfo,2)/bothCam;
        
        for i = 1 : nIt
            idx = 1+(i-1)*bothCam:i*bothCam;
            data = cellfun(@str2double,{frameInfo(idx).T});
            
            [dataSorted,sortIdx] = sort(data);
            
            frameInfo(idx) = frameInfo(sortIdx+idx(1)-1);
            
            if i>1
                if frameInfo(idx(1)).T < frameInfo(idx(1)-1).T
                    
                    dataSorted = dataSorted + (i-1)*(maxFr+1);
                    dataSorted = num2cell(dataSorted);
                    dataSorted = cellfun(@num2str,dataSorted,'UniformOutput',false);
                    for j = 1: length(dataSorted)
                        frameInfo(idx(j)).T = dataSorted{j};
                    end
                    
                    
                end
            end
        end
        for i = 1:size(frameInfo,2)
            if i>1
                currentZ    = str2double(frameInfo(i).Z);
                previousZ   = str2double(frameInfo(i-1).Z);
               

                 if currentZ < previousZ
    % 
                    newZ = currentZ+1;
                    frameInfo(i).T = num2str(newZ);

                 end
                 
                   
%                     currentIFD  = str2double(frameInfo(i).IFD);
%                     if currentIFD ==0
%                     else
%                         previousIFD = str2double(frameInfo(i-1).IFD);
%                         newIFD = previousIFD+1;
%                         frameInfo(i).IFD = num2str(newIFD);
%                     end
            end
%             if i <=2
% 
%             else
%                 currentT    = str2double(frameInfo(i).T);
%                 previousT   = str2double(frameInfo(i-1).T);
%                 previous2T  = str2double(frameInfo(i-2).T);
%                 currentZ    = str2double(frameInfo(i).Z);
%                 previousZ   = str2double(frameInfo(i-1).Z);
%                 currentIFD  = str2double(frameInfo(i).IFD);
%                 previousIFD = str2double(frameInfo(i-1).IFD);
%                 
%                 testDiff = abs(currentT-previousT);
%                 if testDiff>maxFr-2
%                     newT = previousT+1;
%                     frameInfo(i).T = num2str(newT+currentT);
%                 end
%                 
% %                 %test T
% %                 if currentT < previousT
% %                     if previous2T == previousT
% %                         newT = previousT+1;
% % 
% %                     else
% %                         newT = previousT;
% %                     end
% %                     frameInfo(i).T = num2str(newT);
% %                 end
% 
%                

%             end



        end
    end
    assert(totFrame(1) == totFrame(2),'Number of frame of both camera are not consistent');
    %disp('combining done');
end