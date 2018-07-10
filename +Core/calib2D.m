classdef calib2D < Core.Movie
    %This class will hold information about the 2D Cal movie.
    
    properties (SetAccess = 'private')
        cal
    end
    
    methods
        function obj = calib2D(raw)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Core.Movie(raw);
            obj.cal = raw;
          
        end
        
        function set.cal(obj,cal)
            
            if isstruct(cal)
            
            obj.cal = cal;
            
            elseif isfolder(cal)
                
                [file2Analyze] = getFileInPath(obj, cal, '.mat');
                
                if (~isempty(file2Analyze))
                
                    disp('The calibration was already calculated, Loading from existing file');
                    fullPath = [file2Analyze.folder filesep file2Analyze.name];
                    tmp = load(fullPath);
                    cal = tmp.calibration;
                    obj.cal = cal.file;
                    disp('Done');

                    %otherwise we calculate it
                else
                        disp('No calibration Calculating the calibration from calibration data');
                        [calibration] = obj.calcCalibration(cal);
                        obj.cal = calibration;
                        disp('Calibration is now saved');

                end
            else
                error('There is something wrong with your Cal file or path');
           
            end    
        end
        
        function calc(obj)
            
            path = obj.raw.movInfo.Path;
            [file2Analyze] = obj.getFileInPath( path, '.mat');
            
            if (~isempty(file2Analyze))
                
                disp('The calibration was already calculated, Loading from existing file');
                fullPath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullPath);
                calibration = tmp.calibration;

            else
            
                path = obj.raw.fullPath;
                [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(path);
                assert(length(movInfo.Cam)==2,'Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
                assert(length(unique(cellfun(@str2num,{frameInfo.Z})))>2,'Z position is not changing across the selected calibration file, this is strange.');

                [calib, inform] = mpSetup.cali.calculate(fullPath);

                calibration.info = inform;
                calibration.file = calib;

                filename = [file2Analyze.folder filesep 'calibration.mat'];
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

