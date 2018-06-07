classdef Movie <  handle
    %MOVIE will hold most of the information about the data and is the
    %object with which the user will mainly interact.
    
    %   Detailed explanation goes here
    
   properties (SetAccess = 'private')
       
       raw
       cal
       status = 'none';
       info
       calibrated
       superResCal
 
    end
    
    methods
        function obj = Movie(raw, info, cal, calibrated, superResCal)
            %MOVIE Construct an instance of this class
            %   We allow the user to create this object with various number
            %   of input allowing therefore to restart the analysis at any
            %   steps in the process.
            obj.raw = raw;
            %Give a status depending on input.
            switch nargin
                case 1 
                case 2
                    obj.info = info;
                case 3
                    obj.info = info;
                    obj.cal = cal;
                case 4
                     obj.cal = cal;
                     obj.info = info;
                     obj.calibrated = calibrated;
                case 5
                     obj.cal = cal;
                     obj.info = info;
                     obj.calibrated = calibrated;
                     obj.superResCal = superResCal;
                otherwise 
                    error('Unexpected number of argument');
            end
            obj.updateStatus;
        end
        
        function  set.raw(obj,raw)
          
            %Check Given path
            [file2Analyze] = getOMETIF(obj,raw);
      
            if length(file2Analyze)>1
                fprintf('More than one Tiff, Loading %s \n', file2Analyze(1).name);
            end
            fullPath = [file2Analyze.folder filesep file2Analyze(1).name];
            [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(fullPath);
            %Check info for 2 cam
            assert(length(movInfo.Cam)==2,'Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
            obj.raw.movInfo = movInfo;
            obj.raw.frameInfo = frameInfo;
            obj.raw.fullPath = fullPath;
            
            obj.updateStatus;
        end
        
        function set.info(obj,info)
            if checkRes
                obj.info = info;
            else
                warning('Something is wrong with the format of the info provided by the user, storage aborted')
            end
        end
        
        function set.cal(obj,cal)
            assert(obj.checkStatus('raw'),'Experimental status is too early to get the calibration, probably raw is missing');
            [file2Analyze] = getFileInPath(obj, cal, '.mat');
            
            if (~isempty(file2Analyze))
                fullPath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullPath);
                cal = tmp.calibration;
                obj.cal = cal;
            else
                [calibration] = obj.calcCalibration(cal);
                obj.cal = calibration;
            end
            
            
            obj.updateStatus;
        end

        function set.calibrated(obj,calibrated)
          assert(obj.checkStatus('rawCal'),'Experiment status is too early to calibrate, probably calibration and/or raw file missing');
          folderContent = dir(calibrated);
          idx2Calibrated = contains({folderContent.name}, 'calibrated');
          
          if length(unique(idx2Calibrated))<2
              [calibrated] = obj.calibrate;
          elseif length(unique(idx2Calibrated))==2
              fullPath = [calibrated filesep 'calibrated'];
              [file2Analyze] = getFileInPath(obj, fullPath, '.mat');           
              if (~isempty(file2Analyze))
                fullpath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullpath);
                calibrated = tmp.calib;
              else
                [calibrated] = obj.calibrate;
              end
          else
              error('Something is wrong with your calibrated directory');
          end
            obj.calibrated = calibrated;
            obj.updateStatus;
        end
        
        function set.superResCal(obj,superResCal)
            obj.superResCal = superResCal;
        end
        
        function [file2Analyze] = getFileInPath(~, path, ext)
            
            assert(ischar(path),'The given path should be a char');
            assert(ischar(ext),'The given extension should be a char');
            folderContent = dir(path);
            index2Images   = contains({folderContent.name},ext);
            file2Analyze = folderContent(index2Images);
        end
        
        function [calibration] = calcCalibration(obj,path)
            
            [file2Analyze] = obj.getOMETIF(path);
            fullPath = [file2Analyze.folder filesep file2Analyze.name];
            [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(fullPath);
            assert(length(movInfo.Cam)==2,'Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
            assert(length(unique(cellfun(@str2num,{frameInfo.Z})))>2,'Z position is not changing across the selected calibration file, this is strange.');
            
            [calib, inform] = mpSetup.cali.calculate(fullPath);
            
            calibration.info = inform;
            calibration.file = calib;
            
            filename = [file2Analyze.folder filesep 'calibration.mat'];
            calibration.fullPath = filename;
            save(filename,'calibration');
   
        end
        
        function [calib] = calibrate(obj)
           
            frame = 1:obj.raw.movInfo.maxFrame(1);
            [data, ~, ~] = mpSetup.loadAndCal( obj.raw.fullPath, obj.cal.file, frame);
            step = 100;
            calDir = [obj.raw.movInfo.Path filesep 'calibrated'];
            mkdir(calDir);
            
            fid = fopen([calDir filesep 'CalibratedInfo.txt'],'w');
            fprintf(fid,'The information in this file are intended to the user. They are generated automatically so please do not edit them\n');
            calib.mainPath = calDir;
            for i = 1:size(data,3)

                fName = sprintf('calibratedPlane%d.tif',i);
                fPathTiff = [calDir filesep fName];
                fieldN = sprintf('plane%d',i);
                calib.filePath.(fieldN) = fPathTiff;
                t = Tiff(fPathTiff, 'w');
                    for j = 1:step:size(data,4)
                    range = j:j+step-1;
                        if max(range)>= size(data,4)
                        range = j:size(data,4);
                        end
                    t = dataStorage.writeTiff(t,squeeze(data(:,:,i,range)),16);
                    end
                t.close;
                fprintf(fid,...
                'Image plane %d: Cam %d, Channel %d Col1: %d Col2: %d, Rel. Zpos: %0.3f \n ',...
                i,obj.cal.file.inFocus(obj.cal.file.neworder(i)).cam,...
                obj.cal.file.inFocus(obj.cal.file.neworder(i)).ch,...
                obj.cal.file.ROI(obj.cal.file.neworder(i),1),...
                obj.cal.file.ROI(obj.cal.file.neworder(i),1)+...
                obj.cal.file.ROI(obj.cal.file.neworder(i),3),...
                obj.cal.file.inFocus(obj.cal.file.neworder(i)).zpos-...
                obj.cal.file.inFocus(obj.cal.file.neworder(1)).zpos);
                
                
            end
            fclose(fid);
            fName = [calDir filesep 'calibrated.mat'];
            save(fName,'calib');
     
        end
        
        function getCalibration(obj,path)
            obj.cal = path;
        end
        
        function getCalibrated(obj,path)
            obj.calibrated = path;
        end
        
        function [file2Analyze] = getOMETIF(obj,path)
            expExt = '.ome.tif';
            %Check Given path
            [file2Analyze] = obj.getFileInPath(path, expExt);
            assert(~isempty(file2Analyze),sprintf('The given directory does not any %s files',expExt));
        end
        
        function [checkRes] = checkStatus(obj,status)
            %This function check whether the "condition status" is at least
            %reach. e.g. When loading a calibration file we want to make
            %sure that the status is raw or higher.

            checkRes = false;
            switch status
                case 'none'
                        checkRes = true;
                case 'raw'
                    if(strcmp(obj.status,status))
                    checkRes = true;
                    else
                        [checkRes] = obj.checkStatus('rawCal');
                    end
                case 'rawCal'
                    
                    if(strcmp(obj.status,status))
                    checkRes = true;
                    else
                        [checkRes] = obj.checkStatus('calibrated');
                        cond2 = ~isempty(obj.raw);
                        checkRes = logical(checkRes*cond2);
                        
                    end
                    
                case 'calibrated'
                    if(strcmp(obj.status,status))
                        checkRes = true;
                    else    
                        [checkRes] = obj.checkStatus('SRCalibrated');
                    end
                   
                case 'SRCalibrated'
                    if(strcmp(obj.status,status))
                        checkRes = true;
                    end
            end
        end
        
%         function getRaw(obj,raw)
%             if(isempty(obj.raw))
%                 obj.raw = raw;
%             else
%                 obj.raw = raw;
%             end
%             
%         end
%         
%         function getRCal(obj,cal)
%             obj.cal;
%         end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    methods (Access = private)
        
        function updateStatus(obj)
            %Cascade, as long as nothing empty is found we keep looking at
            %the next element to know which status it is. 
            if(isempty(obj.raw))
                obj.status = 'none';
            elseif (isempty(obj.cal))
                obj.status = 'raw';
            elseif (isempty(obj.calibrated))
                obj.status = 'rawCal';
            elseif (isempty(obj.superResCal))
                obj.status = 'calibrated';
            elseif (~isempty(obj.superResCal))
                obj.status = 'SRCalibrated';
            else
                error('Whooop something is wrong with status, go in the corner and think about what you just did!')
            end
            %if the cascade stop in some places we still check that
            %calibrated and SR calibrated are not empty because we allow
            %the user to load from the calibrated data directly(without
            %raw)
            if (~isempty(obj.calibrated))
                if(~isempty(obj.superResCal))
                    obj.status = 'SRCalibrated';
                else
                    obj.status = 'calibrated';
                end
            end

        end
        
    end
end

