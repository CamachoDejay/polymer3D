classdef TrackingExperiment < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        trackMovies
        cal2D
        ZCal
        SRCal
        traces3D
        RMSD
    end
    
    methods
        function obj = TrackingExperiment(folder2Data,cal2D,zCalPath,SRCalPath)
            %TrackingExperiment Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = folder2Data;
            obj.cal2D = cal2D;
            obj.ZCal = zCalPath;
            obj.SRCal = SRCalPath;
        end
        
        function set.path(obj, path)
            assert(ischar(path), 'Path should be given as a string');
            assert(isfolder(path), 'The path given is not a folder, ZCalibration expect a folder. In the folder it is expected to find separate folder for each zCalMovie.')
            
            folderContent = dir(path);
            %Get how many folder are in the main folder
            idx = sum(cellfun(@sum,{folderContent.isdir}));
            %Matlab always store ., .. as folder for relative path so we
            %want to find more than 2 folder in folderContent.
            assert(sum(idx)>2, 'No folder was found in the path given. Expected to find separate folder for each zCalMovie.');
            
            obj.path = path;
            
         end
        
        function set.cal2D(obj,cal2D)
            assert(isstruct(cal2D),'2D Calibration is expected to be received as a structure');
            assert(isfield(cal2D,'fullPath'),'Missing field "fullPath" in cal2D structure');
            assert(isfield(cal2D,'file'),'Missing field "file" in cal2D structure');
            assert(isfield(cal2D,'info'),'Missing field "info" in cal2D structure');
            
            obj.cal2D = cal2D;
            
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
                assert(and(isstruct(cal), and(isfield(cal,'trans'),isfield(cal,'rot'))),...
                    'SR calibration is supposed to be a struct with 2 fields');
                
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
                assert(isstruct(cal),'zCalibration is supposed to be in cells format');
                assert(and(isfield(cal,'fitZParam'),isfield(cal,'calib')),...
                    'Something is wrong in the fields of your Z calibration');
                
                obj.ZCal = cal; 
            end
        end
        
   end
end

