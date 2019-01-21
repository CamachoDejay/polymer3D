classdef MPMovie < Core.Movie
    %From this class we assume that we are dealing with multiplane data.
    %The movie can thus be calibrated if calibration data is provided.
    %Extension will be centered around localization but SOFI movie could
    %also inherit from this object
    
    properties (SetAccess = 'protected')
        cal2D %calibration data
        calibrated %Path to calibrated data and Info
    end
    
    methods
        
        function obj = MPMovie(raw,cal,info)
           
            obj = obj@Core.Movie(raw,info);

            obj.cal2D = cal;
            
        end
        
        function set.cal2D(obj,cal2D)
            
            assert(isstruct(cal2D), 'Calibration is expected to be a struct');
            assert(numel(fieldnames(cal2D))==3, 'Calibration is expected to have 3 Fields');
            assert(isfield(cal2D,'file'),'One of the field should be "file" ');
            
            obj.cal2D = cal2D;
               
        end
        
        function set.calibrated(obj,calibrated)
           
              assert(isstruct(calibrated),'Calibrated is expected to be a struct');
              obj.calibrated = calibrated;
              
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

                   if ~isempty(file) %only work with 8 Planes now
                        %If data is already calibrated we load it
                       disp('The dataset is already calibrated, Loading from existing file');
                       
                       if length(file) ~= 8
                           
                            warning('Did not find 8 planes, if you are not using the prism it is okay, otherwise you might want to recalibrate');
                       
                       end
                   end
                    
                   fullpath = [file2Analyze.folder filesep file2Analyze.name];
                   tmp = load(fullpath);
                   calibrate = tmp.calib;

                    for i = 1: length(fieldnames(calibrate.filePath))

                        currentPath = calibrate.filePath.(['plane' num2str(i)]);
                        idx1 = strfind(fullPath,'\calibrated');
                        idx2 = strfind(currentPath,'\calibrated');                           
                        newPath = [fullPath(1:idx1-1) currentPath(idx2(1):end)];
                        calibrate.filePath.(['plane' num2str(i)]) = newPath;

                    end
                    
                    disp('Done');
                    
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
            [camConfig] = obj.determineCAMConfig;
            obj.calibrated.camConfig = camConfig;

        end
        
        function h = showFrame(obj,idx,scaleBar,idx2Plane)
            switch nargin
                case 3
                    idx2Plane = [];
                case 4
                otherwise
                    error('Not enough input argument');
            end
            
            %Adapted method from the Core.Movie one, its behavior changed
            %depending on whether the data has been calibrated or not.
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            %check the frame requested
            [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
            %Get the data of the requested frame
            [frame] = getFrame(obj,idx);            
            assert(isstruct(frame),'Error unknown data format, data should be a struct');
            
            pxSize = obj.info.pxSize/1000;%in µm
            scaleBarPx = scaleBar/pxSize;
            planes = fields(frame);
            
            if ~isempty(idx2Plane)
                for i = 1 :length(planes)
                    if i ~=idx2Plane
                        frame = rmfield(frame,planes{i});
                    end
                    
                end
                nsFig = 1;
                
            else
                nsFig = 2;
            end
                
            nImages = numel(fields(frame));
            fNames = fieldnames(frame);
            
            
            if ~isempty(obj.calibrated)
                
                zPos = obj.calibrated.oRelZPos;
                
            else
                
                zPos = zeros(size(fNames));
                
            end
            
            %Displaying occur hear
            h = figure(1);
            h.Name = sprintf('Frame %d',idx);
            
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
            
            for i = 1:nImages
                currentIM = frame.(fNames{i});
                subplot(nsFig,nImages/nsFig,i)
                
                if strcmp(obj.info.type,'transmission')
                    colormap('gray');
                    currentIM = imcomplement(currentIM);
                else
                    colormap('jet')
                end
                
                imagesc(currentIM)
                caxis([min(min(min(currentIM))), max(max(max(currentIM)))]);
                hold on 
                x = size(currentIM,2)-scaleBarPx-(0.05*size(currentIM,2)):size(currentIM,2)-0.05*size(currentIM,2);
                y = ones(1,length(x))*size(currentIM,1)-0.05*size(currentIM,2);
                text(mean(x),mean(y)-0.04*size(currentIM,1),[num2str(scaleBar) ' µm'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',8);
                plot(x,y,'-w','LineWidth',2.5);
                axis image;
                
                %grid on;
                a = gca;
                a.XTickLabel = [];
                a.YTickLabel = [];
                a.GridColor = [1 1 1];
                title({fNames{i},sprintf(' Zpos = %0.3f',zPos(i))});
                hold off
               
                
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
        
        function [data] = getPlane(obj,idx,frame)
            %Allow the user to extract data from a specific plane, behavior
            %depends on the calibration
            switch nargin 
                case 1
                    frame = 1:obj.raw.maxFrame(1);
                case 2
            end
            assert(and(idx<9,idx>=1),'plane should be between 1 and 8');
            
            [data] = Load.Movie.tif.getframes(obj.calibrated.filePath.(sprintf('plane%d',idx)),frame);
             
            
        end
        
        function [calibrated] = getCalibrated(obj)
            
            calibrated = obj.calibrated;
            
        end
        
        function playMovie(obj)
            %TODO: code this function
        end
        
        function [xStep, xPosMotor] = getXPosMotor(obj)
              %small function that extract the zStep from the info in the raw
            nFrames = obj.raw.maxFrame(1);
            xPosMotor = zeros(nFrames,1);
            
            for i = 1 : obj.raw.maxFrame(1)
                
                xPosMotor(i) = obj.raw.frameInfo(2*i).Pos(1);
                
            end
            
            xStep = diff(xPosMotor);
            
        end
        
        function [yStep, yPosMotor] = getYPosMotor(obj)
            %small function that extract the zStep from the info in the raw
            nFrames = obj.raw.maxFrame(1);
            yPosMotor = zeros(nFrames,1);

            for i = 1 : obj.raw.maxFrame(1)

                yPosMotor(i) = obj.raw.frameInfo(2*i).Pos(2);

            end

            yStep = diff(yPosMotor);

        end
        
        function [zStep, zPosMotor] = getZPosMotor(obj)
            %small function that extract the zStep from the info in the raw
            nFrames = obj.raw.maxFrame(1);
            zPosMotor = zeros(nFrames,1);
            
            for i = 1 : obj.raw.maxFrame(1)
                
                zPosMotor(i) = obj.raw.frameInfo(2*i).Pos(3);
                
            end
            
            zStep = diff(zPosMotor);
            
        end 
        
        function [camConfig] = determineCAMConfig(obj)
            
             planeDist = abs(mean(diff(obj.calibrated.oRelZPos)))*1000;
             if planeDist > 350
                 camConfig = 'fullRange';
             elseif and(planeDist < 350, planeDist>200)
                 camConfig = 'alternated';
             else
                 error('Something is wrong with your distance between planes')
             end
             
         end
                
    end
    
    methods (Static)
       
    end
    
    methods (Access = private)
        
        function [calib] = applyCalib(obj)
            
            frameInfo = obj.raw.frameInfo;
            movInfo   = obj.raw.movInfo;
            step = 100;
            maxFrame = obj.raw.movInfo.maxFrame(1); 
            frame2Load = obj.info.frame2Load;
            if ischar(frame2Load)
                frame2Load = 1:maxFrame;
            else
                frame2Load = Core.Movie.checkFrame(frame2Load,maxFrame);
                if max(frame2Load) > maxFrame
                   
                    frame2Load = 1:maxFrame;
                                       
                end
            end
            
            endFrame = max(frame2Load);
            startFrame = frame2Load(1);
            
            nStep = ceil((endFrame-startFrame)/step);
            frame = startFrame:step:endFrame;
            nFrame = endFrame-startFrame;
            for i = 1:nStep
                
                if i < nStep
                   
                    cFrame = frame(i):frame(i+1)-1;
                else

                    cFrame = frame(i):endFrame;                    
                    
                end
                 % load the raw data 
                [ movC1, movC2] = Load.Movie.ome.load( frameInfo, movInfo, cFrame );

                %applying the calibration
                [data] = mpSetup.cali.apply( movC1, movC2, obj.cal2D.file );

                %saving data per plane and info to cal
                [calib] = obj.saveCalibrated(data,endFrame);
            end           
         end
         
        function [calib,fid] = saveCalibrated(obj,data,maxFrame)
            cal = obj.cal2D.file;
            calDir = [obj.raw.movInfo.Path filesep 'calibrated'];
            mkdir(calDir);
            %Save the resulting planes in separated TIF and save a txt info
            %file
            fid = fopen([calDir filesep 'CalibratedInfo.txt'],'w');
            fprintf(fid,'The information in this file are intended to the user. They are generated automatically so please do not edit them\n');
            calib.mainPath = calDir;
            calib.nPlanes   = size(data,3);
            
            for i = 1:size(data,3)
                data2Store = squeeze(data(:,:,i,:));
                if strcmpi(obj.info.type,'transmission')
                    data2Store = imcomplement(data2Store);
                end
                
                fName = sprintf('calibratedPlane%d.tif',i);
                fPathTiff = [calDir filesep fName];
                fieldN = sprintf('plane%d',i);
                calib.filePath.(fieldN) = fPathTiff;
                calib.nFrames = maxFrame;
                t = Tiff(fPathTiff, 'a');

                t = dataStorage.writeTiff(t,data2Store,16);

                t.close;
                %We also write a few info about the calibrated data in a
                %text file
                fprintf(fid,...
                'Image plane %d: Cam %d, Channel %d Col1: %d Col2: %d, Rel. Zpos: %0.3f \n ',...
                i,cal.inFocus(cal.neworder(i)).cam,...
                cal.inFocus(cal.neworder(i)).ch,...
                cal.ROI(cal.neworder(i),1),...
                cal.ROI(cal.neworder(i),1)+...
                cal.ROI(cal.neworder(i),3),...
                cal.inFocus(cal.neworder(i)).relZPos);
            
                calib.oRelZPos(i) =  cal.inFocus(cal.neworder(i)).relZPos;
             
            end
            fclose(fid);
            fName = [calDir filesep 'calibrated.mat'];
            save(fName,'calib');
         end

    end
end

