classdef calib2D < Core.Movie
    %This class will hold information about the 2D Cal movie and will be
    %able to calculate the calibration as well.
    
    properties (SetAccess = 'private')
        cal
    end
    
    methods
        function obj = calib2D(raw)
            %In this case raw path is used for both the loading of the
            %movie and the calibration calculation
            obj = obj@Core.Movie(raw);
            obj.calc;
          
        end
        
        function set.cal(obj,cal)
            assert(isstruct(cal),'2DCal expected to be a struct');
            obj.cal = cal;
                
        end
        
        function calc(obj)
            assert(~isempty(obj.raw),'Please provide a path');
            %Calculate from the raw path stored in the movie
            path = obj.raw.movInfo.Path;
            [file2Analyze] = obj.getFileInPath( path, 'calibration.mat');
            
            if (~isempty(file2Analyze))
                %If calibration already exist we load it
                disp('The calibration was already calculated, Loading from existing file');
                fullPath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullPath);
                calibration = tmp.calibration;

            else
                %Otherwise we Calculate it
                path = obj.raw.fullPath;
                [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(path);
                assert(length(movInfo.Cam)==2,'Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
                assert(length(unique(cellfun(@str2num,{frameInfo.Z})))>2,'Z position is not changing across the selected calibration file, this is strange.');

                [calib, inform] = mpSetup.cali.calculate(path);

                calibration.info = inform;
                calibration.file = calib;
                [folder,~] = fileparts(path);
                filename = [folder filesep 'calibration.mat'];
                calibration.fullPath = filename;
                save(filename,'calibration');
                
            end
            
            obj.cal = calibration;
            disp('Done');
            
        end
        
        function cal = getCal(obj)
            
            cal = obj.cal;
            
        end
        
    end
end

