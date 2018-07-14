classdef superResCal < Core.zCalMovie
    %SuperResCal will hold information about the operation that needs to be
    %performed on the coordinates to match them to the plane of reference.
    
    %Although it is not obvious why it should inherit from zCalMovie, the
    %reason is that in essence, the superResCalMovie is a zStack and
    %similar analysis (e.g. localization, consolidation, trackInZ,...) will
    %need to be performed prior to calculate the superResCal.
    
    %Thus rather than having two classes with a lot in common (thus
    %duplicate) it was decided that it would inherit from zCalMovie
    
    properties
        superResCalibration
    end
    
    methods
        %Constructor
        function obj = superResCal(raw, cal, candidatePos)
            
            obj  = obj@Core.zCalMovie(raw,cal);
            
            if nargin == 3
                
                    obj.candidatePos = candidatePos;
                    
            end
        end
        
        
    end
end

