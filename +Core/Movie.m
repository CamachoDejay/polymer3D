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
       
        function set.raw(obj,raw)
            
            assert(isfolder(raw), 'The given path is not a folder');
            %Check Given path
            [file2Analyze] = getOMETIF(obj,raw);
      
            if length(file2Analyze)>1
                
                fprintf('More than one Tiff, Loading %s \n', file2Analyze(1).name);
                
            end
            
            fullPath = [file2Analyze.folder filesep file2Analyze(1).name];
            [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(fullPath);
            %Check info for 2 cam
            assert(length(movInfo.Cam) == 2,'Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
            obj.raw.movInfo   = movInfo;
            obj.raw.frameInfo = frameInfo;
            obj.raw.fullPath  = fullPath;
            obj.raw.maxFrame  = movInfo.maxFrame;
            obj.updateStatus;
            
        end
        
        function set.info(obj,inform)
            
            assert(isstruct(inform),'Information is expected to be a structure');
            names = fieldnames(inform);
          for i = 1:numel(fields(inform))
              
              obj.info.(names{i}) = inform.(names{i});
              
          end
          
          obj.updateStatus;
          
        end
        
        function set.cal(obj,cal)
            
            assert(obj.checkStatus('raw'),'Experimental status is too early to get the calibration, probably raw is missing');
            assert(isfolder(cal), 'The given path is not a folder');
            [file2Analyze] = getFileInPath(obj, cal, '.mat');
            %if there is a .mat in the folder, the calibration was already
            %calculated, we just load it
            
            if (~isempty(file2Analyze))
                
                disp('The calibration was already calculated, Loading from existing file');
                fullPath = [file2Analyze.folder filesep file2Analyze.name];
                tmp = load(fullPath);
                cal = tmp.calibration;
                obj.cal = cal;
                disp('Done');
                
            %otherwise we calculate it
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
          assert(isfolder(calibrated), 'The given path is not a folder');
          folderContent = dir(calibrated);
          idx2Calibrated = contains({folderContent.name}, 'calibrated');
          
          %if there is only 1 diff value, this value must be 0 and thus
          %calibrated folder does not exist thus, we calibrate
          if length(unique(idx2Calibrated))<2
              
              disp('Calibrating the dataset');
              [calibrated] = obj.calibrate;
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
            obj.updateStatus;
            
        end
        
        function set.candidatePos(obj,candidatePos)
            
            assert(iscell(candidatePos), 'Expecting a cell containing candidate positions');
            obj.candidatePos = candidatePos;
            obj.updateStatus;
            
        end
        
        function set.superResCal(obj,superResCal)
            %Assertion are still needed
            obj.superResCal = superResCal;
            obj.updateStatus;
            
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%USER get/set FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%        
       
        function getRaw(obj,path)
            
            obj.raw = path;
            
        end
        
        function getCalibration(obj,path)
            
            obj.cal = path;
            
        end
        
        function calibrate(obj)
            
            assert(obj.checkStatus('rawCal'),'Experiment status is too early to calibrate, probably calibration and/or raw file missing');
            obj.calibrated = obj.raw.movInfo.Path;
            
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
            %Calculate some setup parameters
            sigma_nm = 0.25 * emW/NA;
            FWHMnm = sigma_nm * sqrt(8*log(2));
            FWHM_pix = FWHMnm/pxSize;
            sigmaPix = sigma_nm/pxSize;
            %store info        
            inform.pxSize = pxSize;
            inform.NA = NA;
            inform.emW = emW;
            inform.FWHM_px =  FWHM_pix;
            inform.sigma_px = sigmaPix;
            inform.comment = comment;
            
            obj.info = inform;
            
        end
        
        function findCandidatePos(obj,detectParam, frames)
            
            assert(~isempty(obj.info), 'Missing information about setup to be able to find candidates, please use giveInfo method first');
            
            switch nargin 
                
                case 1
                    
                    run = false;
                    
                case 2
                    
                    frames = 1: obj.calibrated.nFrames;
                    run = true;
                    disp('Running detection on every frame');
                    
                case 3
                    
                    run= true;
                    [frames] = obj.checkFrame(frames);
                    
                otherwise
                    
                    error('too many inputs');
                    
            end
            
            if run
                
                assert(obj.checkStatus('calibrated'),'Data should be calibrated to detect candidate');
                assert(isstruct(detectParam),'Detection parameter should be a struct with two fields');
                nFrames = length(frames);
                currentCandidate = obj.candidatePos;
                
                if(isempty(currentCandidate))
                    
                    candidate = cell(obj.calibrated.nFrames,1);
                    
                else
                    
                    candidate = currentCandidate;
                    
                end
                
                %parameter for localization
                FWHM_pix = obj.info.FWHM_px;
                delta  = detectParam.delta;
                chi2   = detectParam.chi2;
                h = waitbar(0,'detection of candidates...');
                
                for i = 1 : 1:nFrames
                    
                    position = zeros(200,3);
                    [volIm] = obj.getFrame(frames(i)); 
                    nameFields = fieldnames(volIm);
                    
                    for j = 1:length(nameFields)
                        
                        [ pos, ~, ~ ] = Localization.smDetection( double(volIm.(nameFields{j})),...
                            delta, FWHM_pix, chi2 );
                        startIdx = find(position==0,1,'First');
                        pos(:,3) = j;
                        position(startIdx:startIdx+size(pos,1)-1,:) = pos;
                        
                    end
                    
                    idx = find(position==0,1,'First');
                    candidate{frames(i)} = position(1:idx-1,:);
                    waitbar(i/nFrames,h,...
                        sprintf('detection of candidates in Frame %d/%d done',i,nFrames));
                end
                
                close(h);
                obj.candidatePos = candidate;
                
            else
                
                disp('getCandidatePos is a function that detects features in a movie');
                disp('To work, it needs to receive a structure containing 2 detection parameter:');
                disp('delta which is the radius around which detection is performed? usually 6 pixels');
                disp('chi2 which characterize how certain you want to be about the fact that there is a molecule there');
                disp('Usually between 20 and 200');
                
            end
        end
        
        function [candidate] = getCandidatePos(obj, frames)
            
            [idx] = obj.checkFrame(frames);
            candidate = obj.candidatePos{idx};
            
            if isempty(candidate)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
            
            end
        end
        
        function [data] = getFrame(obj,idx)
            
            assert(length(idx)==1,'Requested frame exceed the size of the movie');
            [idx] = obj.checkFrame(idx);
            %Behavior depend on status
            if or(strcmp(obj.status,'raw'),strcmp(obj.status,'rawCalc'))
                %LoadCam
                [movC1,movC2,~] = Load.Movie.ome.load(obj.raw.frameInfo,obj.raw.movInfo,idx);
                data.Cam1 = movC1;
                data.Cam2 = movC2;
                
            elseif or(strcmp(obj.status,'calibrated'),strcmp(obj.status,'SRCalibrated'))
                
                for i = 1:numel(fields(obj.calibrated.filePath))
                    %Load plane
                    [mov] = Load.Movie.tif.getframes(obj.calibrated.filePath.(sprintf('plane%d',i)),idx);
                    data.(sprintf('plane%d',i)) = mov;
                    
                end
            end  
        end
        
        function [particleList] = consolidate(obj,frames,roiSize)
            
            assert(obj.checkStatus('calibrated'),'Data should be calibrated to detect candidate');
            assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');    
            
            switch nargin
                
                case 1

                    frames = 1: obj.calibrated.nFrames;
                    roiSize = 6;
                    disp('Running consolidation on every frame with roi of 6 pixel');

                case 2

                    [frames] = obj.checkFrame(frames);
                    roiSize = 6;
                    disp('no roiSize given, using default radius of 6 pixels')
                    
                case 3
                    
                    [frames] = obj.checkFrame(frames);
                    assert(min(size(roiSize))==1,'RoiSize is expected to be a single number')
                    assert (isnumeric(roiSize),'RoiSize is expected to be a single number');
                    
                otherwise
                    
                    error('Something wrong with input'); 
                    
            end
            
            nFrames = length(frames);
            %allocate for storage
            particleList = cell(1,nFrames);
            
            for i = 1 : 1:nFrames
                
                idx = frames(i);
                [data] = obj.getFrame(idx);
                [candidate] = obj.getCandidatePos(idx);
                
                if isempty(candidate)
                    
                    warning('Frame %d did not contain any candidate',idx);
                    
                else
                    
                    [finalCandidate] = obj.consolidatePos(data, candidate, roiSize);
                    particleList{i} = finalCandidate;
                    
                end
            end          
        end
               
%%%%%%%%%%%%%%%%%%%%%%%%%%% USER DISPLAY FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%

        function showFrame(obj,idx)
            
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            
            [idx] = obj.checkFrame(idx);
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
        
        function showCandidate(obj,idx)
            
            assert(length(idx)==1, 'Only one frame can be displayed at once');
            
            [idx] = obj.checkFrame(idx);
            assert(~isempty(obj.candidatePos{idx}),'There is no candidate found in that frame, check that you ran the detection for that frame');
            
            [frame] = getFrame(obj,idx);
            assert(isstruct(frame),'Error unknown data format, data should be a struct');
            
            nImages = numel(fields(frame));
            fNames = fieldnames(frame);
            nsFig = 2;
            
            candidate = obj.getCandidatePos(idx);
            rowPos    = candidate(:,1);
            colPos    = candidate(:,2);
            planeIdx  = candidate(:,3);
            
            h = figure(2);
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
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET PATH  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [file2Analyze] = getFileInPath(~, path, ext)
           
            assert(ischar(path),'The given path should be a char');
            assert(ischar(ext),'The given extension should be a char');
            assert(isfolder(path),'The path given is not a folder')
            
            folderContent = dir(path);
            index2Images  = contains({folderContent.name},ext);
            file2Analyze  = folderContent(index2Images);
            
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Candidate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [finalCandidate] = consolidatePos(obj, data, frameCandidate, roiSize)
            
            delta = roiSize;
            sig = [obj.info.sigma_px obj.info.sigma_px];
            
            currentk = 1;
            candMet = zeros(length(frameCandidate),7);
            
            for j = 1 : numel(fields(data))
                
                planeCandidate = frameCandidate(frameCandidate(:,3)==j,1:2);
                planeData = data.(sprintf('plane%d',j));

                cM = zeros(length(planeCandidate),6);
                
                for k = 1 : length(planeCandidate)
                    %Get the ROI    
                    [roi_lims] = EmitterSim.getROI(planeCandidate(k,2), planeCandidate(k,1),...
                        delta, size(planeData,2), size(planeData,1));
                    ROI = planeData(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2));
                    %Phasor fitting to get x,y,e
                    [row,col,e] = Localization.phasor(ROI);
                    [LRT,RMSD] = Localization.likelihoodRatioTest(ROI,sig,[row col]);
                    rowPos = planeCandidate(k,1) + row;
                    colPos = planeCandidate(k,2) + col;

                    cM(k,:) = [rowPos colPos e LRT RMSD k];
                
                end
                %candidate Metric: phasor localization Res, Likelihood Res,
                %indexCandidate in plane, Plane.
                %candidateMetric: [row col e LRT RMSD k J]
                
                candMet(currentk:currentk+k-1,:) = [cM j*ones(k,1)];
                currentk = currentk+k;
                
            end
            
            [corrEllip, focusMetric] = Localization.calcFocusMetric(candMet(:,3),candMet(:,4));
            %Kill "bad" PSF from the focus metric so they are not taken
            %later
            focusMetric((1-corrEllip)>0.2) = NaN;
            candMet = [candMet focusMetric];
            
            counter = 1;
            nPart = 0;
            maxIt = size(candMet,1);
            finalCandidate = cell(max(size(find(~isnan(focusMetric)))),2);
            
            while and(~isempty(candMet(:,8)), ~isnan(nanmax(candMet(:,8))))
                
                if counter> maxIt
                    
                    error('While loop ran for an unexpectedly long time, something might be wrong');
               
                end
              
                %Find candidate in best focus
                [~,idx] = max(candMet(:,8));
                
                %Check if there are planes above
                [idx1] = obj.searchPartnerCandidate(candMet,idx,1);

                %Check if they are planes below
                [idx2] = obj.searchPartnerCandidate(candMet,idx,-1);

                idxF = [idx1 idx idx2];
                
                %Check if the resulting configuration of the plane make
                %sense
                [planeConfig] = obj.checkPlaneConfig(idxF,candMet(:,7),idx);
                %Store
                if ~isempty(planeConfig)
                    
                    nPart = nPart +1;
                    finalCandidate{nPart,1} = candMet(idx,:);
                    finalCandidate{nPart,2} = planeConfig;
                    %We remove the particle(s) from the list
                    candMet(idxF,:) = [];
                    
                else
                    %Otherwise we remove it from the best focus search list
                    %by putting focus metric to NaN
                    candMet(idx,8) = NaN;
                    
                end

                counter = counter+1;
                
            end
            %we delete empty cells from the array
            idx2Empty = cellfun(@isempty,finalCandidate);
            finalCandidate(idx2Empty(:,1),:) = [];
            
end
        
        function [idx2Part] = searchPartnerCandidate(obj, candMet, idx, direction)
            
            %This function is designed to have PSFE plate ON
            nP = obj.calibrated.nPlanes; 
            assert(abs(direction) == 1, 'direction is supposed to be either 1 (up) or -1 (down)');
            currentPlane = candMet(idx,7);
            currentPart  = candMet(idx,:);
            
            switch direction

                case 1

                    planes2Check = currentPlane-2:currentPlane-1;
                    planes2Check = planes2Check(planes2Check>0);

                case -1

                    planes2Check = currentPlane+1:currentPlane+2;
                    planes2Check = planes2Check(planes2Check<nP+1);

                otherwise
                    error('direction is supposed to be either 1 (up) or -1 (down)');
            end

            idx2Part = [];
            
            for i = 1 : length(planes2Check)
                
                planeIdx = planes2Check(i);
                partInPlane = candMet(candMet(:,7) == planeIdx,:);
                %Check position %Discuss with Rafa
                partInPlane = partInPlane(and(abs(partInPlane(:,1)-currentPart(1))<2,...
                                              abs(partInPlane(:,2)-currentPart(2)<2)),:);
                if size(partInPlane,1)>1
                    
                    disp('More than one particle in close proximity')
                    
                end
                
                if(isempty(partInPlane))
                    
                    disp('No other candidate in plane with similar position')
                    
                else
                %Check ellipticity 
                    switch direction
                        
                        case 1
                            
                            partInPlane = partInPlane(partInPlane(:,3) > currentPart(3),:);
                            
                        case -1
                            
                            partInPlane = partInPlane(partInPlane(:,3) < currentPart(3),:);
                            
                    end

                    if(isempty(partInPlane))
                        
                        disp('No close candidate in plane consistent with ellipticity')
                        
                    else
                        %check LRT
                        maxExpLRT = abs(currentPart(4))+0.1*abs(currentPart(4));
                        partInPlane = partInPlane(abs(partInPlane(:,4)) < maxExpLRT,:);
                        
                        if(isempty(partInPlane))
                            
                        disp('No close candidate in plane consistent with the LRT parameter')
                        
                        end
                    end
                end

                if ~isempty(partInPlane)
                    
                    newIdx = find(candMet(:,1) == partInPlane(1));
                    
                else
                    
                    newIdx = NaN;

                end
                
                idx2Part = [idx2Part newIdx];
                
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [checkRes] = checkStatus(obj, status)
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
        
        function [frames] = checkFrame(obj, frames) 
            
            testFrame = mod(frames,1);
            
            if all(testFrame<1e-4)
                
            else
                
                frames = round(frames);
                warning('Some frame were not integers, they were rounded');
                
            end
            
            assert(isvector(frames),'Frames should be a vector of integers');
            assert(max(frames) <= obj.raw.maxFrame(1),'Request exceeds max frame');
            assert(min(frames) >0, 'Indexing in matlab start from 1');
            
        end
        
        function [planeConfig] = checkPlaneConfig(~, consIdx, plane, currentIdx)
            %Here we will check that the consolidation found based on the
            %best focused particle make sense with what we would expect and
            %also that we have enough planes.
            assert(length(consIdx) <= size(plane,1),'There is something wrong with your consolidated index and your candidate plane List');
            assert(currentIdx < size(plane,1),'idx to the current candidate exceed the candidate plane list provided');
            
            %Let us test that we have consolidate the particle in at least
            %3 Planes
            test3Planes = length(consIdx) >= 3;
            
            if test3Planes
                
                %we get the index found before and after the currentidx
                %and output a formatted planeConfig
                idx = find(consIdx == currentIdx);
                idxBefore = consIdx(1:idx-1);
                idxBefore = plane(idxBefore);
                
                if length(idxBefore)==1
                    
                    idxBefore = [NaN idxBefore];
                    
                elseif isempty(idxBefore)==1
                    
                    idxBefore = [NaN NaN];
                    
                end
                
                idxAfter  = consIdx(idx+1:end);
                idxAfter = plane(idxAfter);
                
                if length(idxAfter)==1
                    
                    idxAfter = [idxAfter NaN];
                    
                elseif isempty(idxAfter)==1
                    
                    idxAfter = [NaN NaN];
                    
                end
                
                planeConfig = [idxBefore(:)' idxAfter(:)'];
                
            else
                
                planeConfig = [];
                
            end
        end
    end
end
