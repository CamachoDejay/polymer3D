classdef mpMovie < Core.Movie
    %From this class we assume that we are dealing with multiplane data.
    %The movie can thus be calibrated if calibration data is provided.
    %Extension will be centered around localization but SOFI movie could
    %also be a branch of this object
    
    properties (SetAccess = 'private')
        cal2D %calibration data
        calibrated %Path to calibrated data and Info
    end
    
    methods
        
        function obj = mpMovie(raw,cal)
           
            obj = obj@Core.Movie(raw);
            
            switch nargin
                case 1
                case 2
                    obj.cal2D = cal;
                    obj.calibrated = raw;
            end 
        end
        
        function set.cal2D(obj,cal2D)
            
            assert(isstruct(cal2D), 'Calibration is expected to be a struct');
            assert(numel(fieldnames(cal2D))==3, 'Calibration is expected to have 3 Fields');
            assert(isfield(cal2D,'file'),'One of the field should be "file" ');
            
            obj.cal2D = cal2D;
               
        end
        
        function set.calibrated(obj,calibrated)
           
            if ischar(calibrated)

              obj.calibrate;

            else

              assert(isstruct(calibrated),'Calibrated is expected to be a struct');
              obj.calibrated = calibrated;

            end
        end
        
        function calibrate(obj)
            %Method that the user will call to calibrate the data
            folderContent = dir(obj.raw.movInfo.Path);
            idx2Calibrated = contains({folderContent.name}, 'calibrated');

            %if there is only 1 diff value, this value must be 0 and thus
            %calibrated folder does not exist thus, we calibrate
            if length(unique(idx2Calibrated))<2

                disp('Calibrating the dataset');
                [calibrate] = obj.applyCalib;
                disp('Data is now calibrated');

            %if there is 2 value, then calibrated folder exist and then we
            %check if a calibration file is in there.
            elseif length(unique(idx2Calibrated))==2

                fullPath = [obj.raw.movInfo.Path filesep 'calibrated'];
                [file2Analyze] = obj.getFileInPath(fullPath, '.mat'); 

                if (~isempty(file2Analyze))

                    [file] = obj.getFileInPath (fullPath,'.tif');

                    if length(file) == 8 %only work with 8 Planes now
                        %If data is already calibrated we load it
                        disp('The dataset is already calibrated, Loading from existing file');
                        fullpath = [file2Analyze.folder filesep file2Analyze.name];
                        tmp = load(fullpath);
                        calibrate = tmp.calib;
                        disp('Done');

                    else
                    %Otherwise we apply the calibration to the data
                    disp('Some planes are missing (expect 8), recalibrating...');
                    [calibrate] = obj.applyCalib;
                    disp('Data is now correctly calibrated');

                    end

                else
                %If no calibration was found we calibrate again
                disp('Calibrating the dataset');
                [calibrate] = obj.applyCalib;
                disp('Data is now calibrated');

                end
            else

              error('Something is wrong with your calibrated directory');

            end

            obj.calibrated = calibrate;

        end
        
        function showFrame(obj,idx)
            %Adapted method from the Core.Movie one, its behavior changed
            %depending on whether the data has been calibrated or not.
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            %check the frame requested
            [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
            %Get the data of the requested frame
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
            %Displaying occur hear
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
            %Allow the user to extract data from a specific frame, behavior
            %depends on the calibration
            assert(length(idx)==1,'Requested frame exceed the size of the movie');
             [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
            %Behavior depend on status
            if isempty(obj.calibrated)
                
               [data] = getFrame@Core.Movie(obj,idx);
                
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
        
        function playMovie(obj)
            %TODO: code this function
        end
        
        function [zStep, zPosMotor] = getZPosMotor(obj)
            %small function that extract the zStep from the info in the raw
            nFrames = obj.raw.maxFrame(1);
            zPosMotor = zeros(nFrames,1);
            
            for i = 1 : obj.raw.maxFrame(1)
                
                zPosMotor(i) = obj.raw.frameInfo(2*i).Pos(3);
                
            end
            
            zStep = mean(diff(zPosMotor));
            
        end
                
    end
    
    methods (Access = private)
        
         function [calib] = applyCalib(obj)
            %Load and calibrate the movie using the calibration file
            frame = 1:obj.raw.movInfo.maxFrame(1);
            [data, ~, ~] = mpSetup.loadAndCal( obj.raw.fullPath, obj.cal2D.file, frame);
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
                    %Saving occurs in 100 steps to avoid memory issue in
                    %matlab
                    for j = 1:step:calib.nFrames
                        
                    range = j:j+step-1;
                    
                        if max(range)>= calib.nFrames
                            
                        range = j:calib.nFrames;
                        
                        end
                        
                    t = dataStorage.writeTiff(t,squeeze(data(:,:,i,range)),16);
                    
                    end
                    
                t.close;
                %We also write a few info about the calibrated data in a
                %text file
                fprintf(fid,...
                'Image plane %d: Cam %d, Channel %d Col1: %d Col2: %d, Rel. Zpos: %0.3f \n ',...
                i,obj.cal2D.file.inFocus(obj.cal2D.file.neworder(i)).cam,...
                obj.cal2D.file.inFocus(obj.cal2D.file.neworder(i)).ch,...
                obj.cal2D.file.ROI(obj.cal2D.file.neworder(i),1),...
                obj.cal2D.file.ROI(obj.cal2D.file.neworder(i),1)+...
                obj.cal2D.file.ROI(obj.cal2D.file.neworder(i),3),...
                obj.cal2D.file.inFocus(obj.cal2D.file.neworder(i)).zpos-...
                obj.cal2D.file.inFocus(obj.cal2D.file.neworder(1)).zpos);
                calib.oRelZPos(i) =  obj.cal2D.file.inFocus(obj.cal2D.file.neworder(i)).zpos-...
                obj.cal2D.file.inFocus(obj.cal2D.file.neworder(1)).zpos;
                
            end
            %Finally we save the calibratedInfo as a matfile
            fclose(fid);
            fName = [calDir filesep 'calibrated.mat'];
            save(fName,'calib');
            
         end

    end
end
