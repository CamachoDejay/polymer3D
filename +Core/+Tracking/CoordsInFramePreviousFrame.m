classdef CoordsInFramePreviousFrame
   properties
       data  
       time
   end
   methods
       function obj=CoordsInFramePreviousFrame()
          obj.data =[];
          obj.time =[];
           
       end
   %% Constructs memory array for next tracking run
       function ArrayToBeginTracking =AddMemoryToPreviouslyTrackedData(obj,LastFrame,MemoryArray_data)
               ArrayToBeginTracking = [LastFrame ; MemoryArray_data];
       end
   end
end