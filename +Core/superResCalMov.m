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
        function obj = superResCalMov(raw, cal, candidatePos)
            
            obj  = obj@Core.zCalMovie(raw,cal);
            
            if nargin == 3
                
                    obj.candidatePos = candidatePos;
                    
            end
        end
        
        function superResCalibrate(obj)
            
            %#1 Track particle in Z
            [traces, counter] = obj.trackInZ;
            
            %#2 Extract Data per particles
            [partData] = obj.extractPartData;
            %#2 Find frames where a particle is approx equally defocused in
            %2 planes.
            [test]= obj.findDefocusedFrame(traces);
        end
    end
    
    methods (Access = private)
       
        function [partData] = extractPartData(obj)
            list = obj.particles.List;
            traces = obj.particles.Traces;
            nTraces = obj.particles.nTraces;
            
            %Allocate Memory
            
            partData = cell(1,nTraces);
            
            for i = 1 : length(traces)
                if ~isempty(traces{i})
                    for j = 1 : length(traces{i})
                        idx2Part = traces{i}{j};
                        if ~isnan(idx2Part)
                        startIdx = length(partData{idx2Part})+1;
                        partData{idx2Part}(startIdx,:) = [list{i}{j}(3,[1:3, end]) i];
                        else
                        end
                    end
                end
            end         
        end
        
       function [defocusFramePer] = findDefocusedFrame(traces)
           
           
       end
        
    end
    
    
end

