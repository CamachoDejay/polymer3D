classdef CoordsInFrameNextFrame
    properties
        dataNext
        timeNext 
        PossibleNeighbours
    end
    
    methods 
        function obj = CoordsInFrameNextFrame()
            obj.dataNext =[];
            obj.timeNext =[];
            obj.PossibleNeighbours=[];
            
            
            
        end
    %% Constructs a cell array out of all the particle ids from the previous image
        function AllIds = ConstructACellArrayWithAllParticleIds(obj,IMPREV,AllIds)
            if ~isempty(IMPREV)
                    AllIds{end+1} = {IMPREV(1,end) };
                    IMPREV(1,:) =[];
                    AllIds = ConstructACellArrayWithAllParticleIds(obj,IMPREV,AllIds);
            end
        end
        
    end
    
    
end