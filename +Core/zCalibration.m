classdef zCalibration < handle
    %zCalibration is a class that will hold the information of several
    %zCalMovie and be able to process these information to create a
    %calibration curve together with displaying etc...
    
    properties
        path
        zCalMovies
        calib
        
    end
    
    methods
        function obj = zCalibration(pathzCal,pathCal2D)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path.zCalMovies = pathzCal;
            obj.path.cal2D = pathCal2D;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

