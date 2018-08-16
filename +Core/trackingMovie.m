classdef TrackingMovie < Core.MPLocMovie
    %trackingMovie Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        traces
        zCalibration
    end
    
    methods
        function obj = TrackingMovie(raw, cal, zCalibration)
            %trackingMovie Construct an instance of this class
            %   Detailed explanation goes here
            obj  = obj@Core.MPLocMovie(raw,cal);
            
            obj.zCalibration = zCalibration;
        end
        
        function set.zCalibration(obj, zCalibration)
            assert(isstruct(zCalibration),'zCalibration is expected to be a struct with 4 fields')
            assert(and(and(isfield(zCalibration,'file'),isfield(zCalibration,'data')),...
                and(isfield(zCalibration,'fitZParam'),isfield(zCalibration,'trackParam'))),...
                'zCalibration is expected to be a struct with 4 fields');
            
            obj.zCalibration = zCalibration;
            
        end
        
        function giveZCal(obj,zCalibration)
            
            obj.zCalibration = zCalibration;
            
        end
                
    end
end

