classdef superResCalMov < Core.zCalMovie
    %SuperResCal will hold information about the operation that needs to be
    %performed on the coordinates to match them to the plane of reference.
    
    %Although it is not obvious why it should inherit from zCalMovie, the
    %reason is that in essence, the superResCalMovie is a zStack and
    %similar analysis (e.g. localization, consolidation, trackInZ,...) will
    %need to be performed prior to calculate the superResCal.
    
    %Thus rather than having two classes with a lot in common (thus
    %duplicate) it was decided that it would inherit from zCalMovie
    
    properties (SetAccess = 'private')
        superResCalibration
    end
    
    methods
        %Constructor
        function obj = superResCalMov(raw, cal)
            
            obj  = obj@Core.zCalMovie(raw,cal);
                        
        end
        
        function superResCalibrate(obj,trackParam)
            
            %#1 Track particle in Z (= Consolidation between frames)
            [traces, counter] = obj.trackInZ(trackParam);
            
            %#2 Extract Data per particles
            [partData] = obj.extractPartData;
            
            %#2 Find frames where a particle is approx equally defocused in
            %2 planes.
            [test]= obj.findDefocusedFrame(partData);
            
        end
    end
    
    methods (Access = private)
       
        function [partData] = extractPartData(obj)
            list = obj.particles.List;
            traces = obj.traces.trace;
            nTraces = obj.traces.nTrace;
            
            %Allocate Memory
            partData = cell(1,nTraces);
            %Extraction of particule data occurs here
            for i = 1 : length(traces)
                if ~isempty(traces{i})
                    for j = 1 : length(traces{i})
                        idx2Part = traces{i}{j};
                        if ~isnan(idx2Part)
                        startIdx = size(partData{idx2Part},1)+1;
                        currentData = list{i}{j}(~isnan(list{i}{j}(:,1)),[1:3, end]);
                        idx = repmat(i,size(currentData,1),1);
                        partData{idx2Part}(startIdx:startIdx+size(currentData,1)-1,:) = [list{i}{j}(~isnan(list{i}{j}(:,1)),[1:3, end]) idx];
                        else
                        end
                    end
                end
            end
            
            % Here we clean the data keeping only the particles that were
            % observed in every planes
            for i = 1: size(partData,2)
                
                data2Test = partData{i};
                test = unique(data2Test(:,4));
                
                if length(test)<obj.calibrated.nPlanes
                    %put empty cells where test fails
                    partData{i} = [];
                end
                
                
            end
            %Deleting empty cells of the cell array
            partData(cellfun('isempty',partData)) = []; 
           
        end
        
        function [defocusFrame] = findDefocusedFrame(obj,partData)
           
            defocusFrame = cell(size(partData));
            
            for i = 1:size(partData,2)
                currentData = partData{i};
                planes = unique(currentData(:,4));
                for j = 1:max(planes)-1
                    dataPlaneA = currentData(currentData(:,4) == planes(j),:);
                    dataPlaneB = currentData(currentData(:,4) == planes(j+1),:);
                    
                    idx = obj.findOptimalDefocusing(dataPlaneA,dataPlaneB);
                    
                    defocusFrame{i}(j,:) = [j idx]; 
                    
                    
                end
            end
            
           
        end
       
        function [idx] = findOptimalDefocusing(obj,dataPlaneA,dataPlaneB)
            nFrames = unique(dataPlaneA(dataPlaneA(:,3)>1,end));
            bestVal = [2 1];
            for i = 1: length(nFrames)
                currentFrame = nFrames(i);
                currentVal = [dataPlaneA(dataPlaneA(:,5)==currentFrame,3) dataPlaneB(dataPlaneB(:,5)==currentFrame,3)];
                cmpVal = [abs(1-1/currentVal(1)) abs(1-currentVal(2))];
                
                if  and(abs(cmpVal(1) - cmpVal(2))< abs(bestVal(1) - bestVal(2)), currentVal(2)<1)
                    bestVal = cmpVal;
                    idx = nFrames(i);
                end
                
            end
            
        end
        
    end
    
    
end

