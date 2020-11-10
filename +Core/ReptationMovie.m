classdef ReptationMovie < Core.MPMovie
    %ReptationMovie a class to analyze movie of molecules reptating 
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
        
        function [data] = cropAllFrames(obj,ROI)
            ROI = round(ROI);
            nFrames = obj.raw.movInfo.maxFrame(1);
            
            h = waitbar(0,'Cropping data');
            data = zeros(ROI(4),ROI(3),obj.calibrated.nPlanes,nFrames);
            for i = 1:nFrames
               frame = obj.getFrame(i);
               
               data(:,:,:,i) = frame(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1,:); 
               
               waitbar(i/nFrames,h,'Cropping data');
            end
            close(h)
            
        end
      
    end
end

