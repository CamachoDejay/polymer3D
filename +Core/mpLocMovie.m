classdef MPLocMovie < Core.MPParticleMovie
    %MPLocMovie hold all the method and info necessary for localization,
    %expension of this object will be first centered around Tracking but
    %it will also be easy to extend to STORM/PALM.
    
    properties (SetAccess = 'private')
        
        SRCal
        ZCal
        corrLocPos
        
    end
    
    methods
        
        function obj = MPLocMovie(raw, MPCal, SRCal, zCal)
            
            obj  = obj@Core.MPParticleMovie(raw,MPCal);
            
            switch nargin
                
                case 1
                    error('MPCal is required to create MPLocMovie')
                case 2
                    obj.SRCal = [];
                    obj.ZCal  = [];
                case 3
                    obj.SRCal = SRCal;
                    obj.ZCal = [];
                case 4
                    obj.SRCal = SRCal;
                    obj.ZCal = zCal;
                otherwise
                    error('Too many input arguments');
            end
        end
        
        function set.SRCal(obj,SRCal)
            assert(isfolder(SRCal), 'The given path is not a folder');
            
            %Check Given path
            [file2Analyze] = Core.Movie.getFileInPath(SRCal,'SRCalibration.mat');
            
            if isempty(file2Analyze)
                error('No SR calibration file found in the given folder');
            else
                fileName = [file2Analyze.folder filesep file2Analyze.name];
                cal = load(fileName);
                field = fieldnames(cal);
                cal = cal.(field{1});
                assert(and(isfield(cal,'translation'),isfield(cal,'rotation')),...
                    'SR calibration is supposed to have a tranlsation and rotation field');
                
                obj.SRCal = cal; 
            end
            
        end
        
        function set.ZCal(obj,zCal)
            
            assert(isfolder(zCal), 'The given path is not a folder');
            
            %Check Given path
            [file2Analyze] = Core.Movie.getFileInPath(zCal,'zCalibration.mat');
            
            if isempty(file2Analyze)
                error('No z calibration file found in the given folder');
            else
                fileName = [file2Analyze.folder filesep file2Analyze.name];
                cal = load(fileName);
                field = fieldnames(cal);
                cal = cal.(field{1});
                assert(iscell(cal),'zCalibration is supposed to be in cells format');
                assert(length(cal) == obj.calibrated.nPlanes,'Unexpected number of fit in zCalibration');
                
                obj.ZCal = cal; 
            end
        end
        
        
    end
    
    methods (Access = private)
  
    end
end