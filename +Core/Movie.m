classdef Movie < handle
    %MOVIE will hold most of the information about the data and is the
    %object with which the user will mainly interact.
    
    %   Detailed explanation goes here
    
   properties (SetAccess = 'private')
       
       raw
       cal
       status
       info
       calibrated
       superResCal
 
    end
    
    methods
        function obj = Movie(raw, cal, info, calibrated, superResCal)
            %MOVIE Construct an instance of this class
            %   Detailed explanation goes here
            assert(and(isstruct(raw),numel(fieldnames(raw))<3),...
                'raw should be a structure with no more than 2 fields');
            
            assert(and(isstruct(cal),numel(fieldnames(cal))<4),...
                'Cal should be a structure with no more than 3 fields'); 
            
            obj.raw = raw;
            obj.cal = cal;
      
            %Give a status depending on input.
            switch nargin
                case 1
                    obj.status = 'raw';
                case 2
                    obj.status = 'rawCal';
                    
                case 3
                    obj.status = 'rawCal';
                    obj.info = info;
                case 4
                     obj.status = 'rawCal';
                     obj.info = info;
                     obj.calibrated = calibrated;
                case 5
                     obj.status = 'calibrated';
                     obj.superResCal = superResCal;
                     obj.info = info;
                     obj.calibrated = calibrated;
                otherwise 
                    error('Unexpected number of argument');
            end
        end
                
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

