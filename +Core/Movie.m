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
       candidatePos
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
%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        function set.info(obj,inform)
            names = fieldnames(inform);
          for i = 1:numel(fields(inform))
              obj.info.(names{i}) = inform.(names{i});
          end
          
        end
        
        function set.cal(obj,cal)
            assert(obj.checkStatus('raw'),'Experimental status is too early to get the calibration, probably raw is missing');
            [file2Analyze] = getFileInPath(obj, cal, '.mat');
            
            if (~isempty(file2Analyze))
                disp('The calibration was already calculated, Loading from existing file');
                fullPath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullPath);
                cal = tmp.calibration;
                obj.cal = cal;
                disp('Done');
            else
                disp('Calculating the calibration from calibration data');
                [calibration] = obj.calcCalibration(cal);
                obj.cal = calibration;
                disp('Calibration is now saved');
            end
            
            
            obj.updateStatus;
        end

        function set.calibrated(obj,calibrated)
          assert(obj.checkStatus('rawCal'),'Experiment status is too early to calibrate, probably calibration and/or raw file missing');
          folderContent = dir(calibrated);
          idx2Calibrated = contains({folderContent.name}, 'calibrated');
          
          if length(unique(idx2Calibrated))<2
              disp('Calibrating the dataset');
              [calibrated] = obj.calibrate;
              disp('Data is now calibrated');
          elseif length(unique(idx2Calibrated))==2
              fullPath = [calibrated filesep 'calibrated'];
              [file2Analyze] = getFileInPath(obj, fullPath, '.mat');           
              if (~isempty(file2Analyze))
                disp('The dataset is already calibrated, Loading from existing file');
                fullpath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullpath);
                calibrated = tmp.calib;
                disp('Done');
              else
                disp('Calibrating the dataset');
                [calibrated] = obj.calibrate;
                disp('Data is now calibrated');
              end
          else
              error('Something is wrong with your calibrated directory');
          end
            obj.calibrated = calibrated;
            obj.updateStatus;
        end
        
        function set.candidatePos(obj,candidatePos)
            assert(iscell(candidatePos), 'Expecting a cell containing candidate positions');
            
            obj.candidatePos = candidatePos;
        end
        
        function set.superResCal(obj,superResCal)
            obj.superResCal = superResCal;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET PATH  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [file2Analyze] = getFileInPath(~, path, ext)
            
            assert(ischar(path),'The given path should be a char');
            assert(ischar(ext),'The given extension should be a char');
            folderContent = dir(path);
            index2Images   = contains({folderContent.name},ext);
            file2Analyze = folderContent(index2Images);
        end
        
        function [file2Analyze] = getOMETIF(obj,path)
            expExt = '.ome.tif';
            %Check Given path
            [file2Analyze] = obj.getFileInPath(path, expExt);
            assert(~isempty(file2Analyze),sprintf('The given directory does not any %s files',expExt));
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calibration  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GETTER FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function getRaw(obj,path)
            obj.raw = path;
        end
        
        function getCalibration(obj,path)
            obj.cal = path;
        end
        
        function getCalibrated(obj,path)
            obj.calibrated = path;
        end
        
        function giveInfo(obj)
            prompt = {'Enter the pixel size: ','Enter the NA of the objective ',...
                'Enter the emission wavelength', 'Any comment about experiment?'};
            dlgTitle = 'Information about experimental parameters';
            numLines = 1;
            defaultVal = {'95','1.2','520',''};
            answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);
            
            assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')

            pxSize = str2double(answer(1));
            assert(~isnan(pxSize),'Number of Frame should be numerical');%If not a number
            
            NA = str2double(answer(2));
            assert(~isnan(NA),'NA should be numerical');
            
            emW = str2double(answer(3));
            assert(~isnan(emW),'Emission wavelength should be numerical');
            
            comment = answer(4);
            
            inform.pxSize = pxSize;
            inform.NA = NA;
            inform.emW = emW;
            inform.comment = comment;
            obj.info = inform;
        end
        
        function getCandidatePos(obj,detectParam, frames)
            switch nargin
                case 3
                    assert(obj.checkStatus('calibrated'),'Data should be calibrated to detect candidate');
                    assert(isvector(frames),'Frames should be a vector of integers');
                    assert(isstruct(detectParam),'Detection parameter should be a struct with two fields');
                    nFrames = length(frames);

                    currentCandidate = obj.candidatePos;
                    if(isempty(currentCandidate))
                        candidate = cell(obj.calibrated.nFrames,1);
                    else
                        candidate = currentCandidate;
                    end

                    %parameter for localization
                    objNA  = obj.info.NA;
                    emWave = obj.info.emW;
                    sigma_nm = 0.25 * emWave/objNA;
                    FWHMnm = sigma_nm * sqrt(8*log(2));
                    pxSnm  = obj.info.pxSize;
                    FWHM_pix = FWHMnm/pxSnm;
                    delta  = detectParam.delta;
                    chi2   = detectParam.chi2;

                    position = zeros(200,3);
                    for i = 1 : 1:nFrames
                        [volIm] = obj.getFrame(frames(i)); 
                        nameFields = fieldnames(volIm);
                        for j = 1:length(nameFields)
                            [ pos, ~, ~ ] = Localization.smDetection( double(volIm.(nameFields{j})),...
                                delta, FWHM_pix, chi2 );
                            startIdx = find(position==0,1,'First');
                            pos(:,3) = j;
                            position(startIdx:startIdx+size(pos,1)-1,:) = pos;
                        end
                        candidate{frames(i)} = position;
                    end

                    obj.candidatePos = candidate;
                otherwise
                    disp('getCandidatePos is a function that detects features in a movie');
                    disp('To work, it needs to receive a structure containing 2 detection parameter:');
                    disp('delta which is the radius around which detection is performed? usually 6 pixels');
                    disp('chi2 which characterize how certain you want to be about the fact that there is a molecule there');
                    disp('Usually between 20 and 200');
            end
        end
        
        function [data] = getFrame(obj,idx)
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            if or(strcmp(obj.status,'raw'),strcmp(obj.status,'rawCalc'))
                [movC1,movC2,~] = Load.Movie.ome.load(obj.raw.frameInfo,obj.raw.movInfo,idx);
                data.Cam1 = movC1;
                data.Cam2 = movC2;
                
            elseif or(strcmp(obj.status,'calibrated'),strcmp(obj.status,'SRCalibrated'))
                for i = 1:numel(fields(obj.calibrated.filePath))
                    [mov] = Load.Movie.tif.getframes(obj.calibrated.filePath.(sprintf('plane%d',i)),idx);
                    data.(sprintf('plane%d',i)) = mov;
                end
            end
                
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%% DISPLAY FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function showFrame(obj,idx)
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            [frame] = getFrame(obj,idx);
            assert(isstruct(frame),'Error unknown data format, data should be a struct');
            nImages = numel(fields(frame));
            fNames = fieldnames(frame);
            nsFig = 2;
            if obj.checkStatus('calibrated')
                zPos = obj.calibrated.oRelZPos;
            else
                zPos = zeros(size(fNames));
            end
            h = figure(10);
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
        
        function showCandidate(obj,idx)
            assert(length(idx)==1, 'Only one frame can be displayed at once');
            assert(~isempty(obj.candidatePos{idx}),'There is no candidate found in that frame, check that you ran the detection for that frame');
            
            [frame] = getFrame(obj,idx);
            assert(isstruct(frame),'Error unknown data format, data should be a struct');
            
            nImages = numel(fields(frame));
            fNames = fieldnames(frame);
            nsFig = 2;
            rowPos = obj.candidatePos{idx}(:,1);
            colPos = obj.candidatePos{idx}(:,2);
            planeIdx = obj.candidatePos{idx}(:,3);
            h = figure(10);
            h.Name = sprintf('Frame %d',idx);
            for i = 1:nImages
                
                subplot(2,nImages/nsFig,i)
                hold on
                imagesc(frame.(fNames{i}))
                scatter(colPos(planeIdx==i),rowPos(planeIdx==i))
                axis image;
                grid on;
                a = gca;
                a.XTickLabel = [];
                a.YTickLabel = [];
                a.GridColor = [1 1 1];
                title({fNames{i},sprintf(' Zpos = %0.3f',obj.calibrated.oRelZPos(i))});
                colormap('jet')
                hold off
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%                
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

