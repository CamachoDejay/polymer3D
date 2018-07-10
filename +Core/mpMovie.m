classdef mpMovie < Core.Movie
    %This class will hold the path to the movie and the calibrated movie.
    
    properties
        cal2D
        calibrated
    end
    
    methods
        
        function obj = mpMovie(raw,cal)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Core.Movie(raw);
            
            switch nargin
                case 1
                case 2
                    obj.cal = cal;
            end 
        end
        
        function set.cal2D(obj,cal2D)
            
            assert(isstruct(cal2D), 'Calibration is expected to be a struct');
            assert(numel(fieldnames(cal2D))==3, 'Calibration is expected to have 3 Fields');
            assert(isfield(cal2D,'file'),'One of the field should be "file" ');
            
            obj.cal = cal2D;
               
        end
        
        function set.calibrated(obj,calibrated)
            
          assert(isfolder(calibrated), 'The given path is not a folder');
          folderContent = dir(calibrated);
          idx2Calibrated = contains({folderContent.name}, 'calibrated');
          
          %if there is only 1 diff value, this value must be 0 and thus
          %calibrated folder does not exist thus, we calibrate
          if length(unique(idx2Calibrated))<2
              
              disp('Calibrating the dataset');
              [calibrated] = obj.applyCalib;
              disp('Data is now calibrated');
              
          %if there is 2 value, then calibrated folder exist and then we
          %check if a calibration file is in there.
          elseif length(unique(idx2Calibrated))==2
              
              fullPath = [calibrated filesep 'calibrated'];
              [file2Analyze] = getFileInPath(obj, fullPath, '.mat'); 
              
              if (~isempty(file2Analyze))
                  
                [file] = getFileInPath (obj,fullPath,'.tif');
                
                if length(file) == 8
                    
                disp('The dataset is already calibrated, Loading from existing file');
                fullpath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullpath);
                calibrated = tmp.calib;
                disp('Done');
                
                else
                    
                    %error('Some planes are missing (expect 8), recalibrating...');
                    disp('Some planes are missing (expect 8), recalibrating...');
                    [calibrated] = obj.applyCalib;
                    disp('Data is now calibrated');
                    
                end
                
              else
                  
                disp('Calibrating the dataset');
                [calibrated] = obj.applyCalib;
                disp('Data is now calibrated');
                
              end
          else
              
              error('Something is wrong with your calibrated directory');
              
          end
          
            obj.calibrated = calibrated;
            
        end
        
        function calibrate(obj)
            
            obj.calibrated = obj.raw.movInfo.Path;
            
        end
        
        function showFrame(obj,idx)
            
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            
            [idx] = obj.checkFrame(idx);
            [frame] = getFrame(obj,idx);            
            assert(isstruct(frame),'Error unknown data format, data should be a struct');
            
            nImages = numel(fields(frame));
            fNames = fieldnames(frame);
            nsFig = 2;
            
            if ~isempty(obj.calibrated)
                
                zPos = obj.calibrated.oRelZPos;
                
            else
                
                zPos = zeros(size(fNames));
                
            end
            h = figure(1);
            h.Name = sprintf('Frame %d',idx);
            for i = 1:nImages
                
                subplot(2,nImages/nsFig,i)
                imagesc(frame.(fNames{i}))
                axis image;
                grid on;
                a = gca;
                a.XTickLabel = [];
                a.YTickLabel = [];
                a.GridColor = [1 1 1];
                title({fNames{i},sprintf(' Zpos = %0.3f',zPos(i))});
                colormap('jet')
                
            end
        end
        
        function [data] = getFrame(obj,idx)
            
            assert(length(idx)==1,'Requested frame exceed the size of the movie');
            [idx] = obj.checkFrame(idx);
            %Behavior depend on status
            if isempty(obj.calibrated)
                %LoadCam
                [movC1,movC2,~] = Load.Movie.ome.load(obj.raw.frameInfo,obj.raw.movInfo,idx);
                data.Cam1 = movC1;
                data.Cam2 = movC2;
                
            elseif isstruct(obj.calibrated)
                
                for i = 1:numel(fields(obj.calibrated.filePath))
                    %Load plane
                    [mov] = Load.Movie.tif.getframes(obj.calibrated.filePath.(sprintf('plane%d',i)),idx);
                    data.(sprintf('plane%d',i)) = mov;
                    
                end
            end  
        end
        
        function [calibrated] = getCalibrated(obj)
            
            calibrated = obj.calibrated;
            
        end
                
    end
    
    methods (Access = private)
        
         function [calib] = applyCalib(obj)
            %Load and calibrate the movie using the calibration file
            frame = 1:obj.raw.movInfo.maxFrame(1);
            [data, ~, ~] = mpSetup.loadAndCal( obj.raw.fullPath, obj.cal.file, frame);
            step = 100;
            calDir = [obj.raw.movInfo.Path filesep 'calibrated'];
            mkdir(calDir);
            %Save the resulting planes in separated TIF and save a txt info
            %file
            fid = fopen([calDir filesep 'CalibratedInfo.txt'],'w');
            fprintf(fid,'The information in this file are intended to the user. They are generated automatically so please do not edit them\n');
            calib.mainPath = calDir;
            calib.nPlanes   = size(data,3);
            for i = 1:size(data,3)

                fName = sprintf('calibratedPlane%d.tif',i);
                fPathTiff = [calDir filesep fName];
                fieldN = sprintf('plane%d',i);
                calib.filePath.(fieldN) = fPathTiff;
                calib.nFrames = size(data,4);
                t = Tiff(fPathTiff, 'w');
                
                    for j = 1:step:calib.nFrames
                        
                    range = j:j+step-1;
                    
                        if max(range)>= calib.nFrames
                            
                        range = j:calib.nFrames;
                        
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
                calib.oRelZPos(i) =  obj.cal.file.inFocus(obj.cal.file.neworder(i)).zpos-...
                obj.cal.file.inFocus(obj.cal.file.neworder(1)).zpos;
                
            end
            
            fclose(fid);
            fName = [calDir filesep 'calibrated.mat'];
            save(fName,'calib');
            
         end

    end
end

