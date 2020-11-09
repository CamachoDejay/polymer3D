classdef ReptationMovie < Core.MPMovie
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = ReptationMovie(raw,cal,info)
            %ReptationMovie Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Core.MPMovie(raw,cal,info);
        end
        
      
    end
end

