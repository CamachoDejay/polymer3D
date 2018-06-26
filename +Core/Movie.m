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
       particles
       
       superResCal
       
   end
   
    properties (SetAccess = 'public')
    zCalibration
       
   
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
        
        function set.particles(obj,particles)
            obj.particles = particles;
        end
        
        function set.zCalibration(obj, zCalibration)
            obj.zCalibration = zCalibration;
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
            [file2Analyze] = getFileInPath(obj, obj.raw.movInfo.Path, '.mat');
            
            if any(contains({file2Analyze.name},'candidatePos')==true)
                quest = 'Some candidate were found in the raw folder, do you want to load them or run again ?';
                title = 'Question to User';
                btn1  = 'Load';
                btn2 = 'run again';
                defbtn = 'Load';
                answer = questdlg(quest,title,btn1,btn2,defbtn);
                switch answer
                case 'Load'
                    candidate = load([file2Analyze(1).folder filesep 'candidatePos.mat']);
                    candidate = candidate.candidate;
                    run = false;
                case 'run again'
                    run = true;
                end
            else
                run = true;
            end    
            
            if run
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

                        fileName = sprintf('%s%scandidatePos.mat',obj.raw.movInfo.Path,'\');
                        save(fileName,'candidate');
                        


                    else

                        disp('getCandidatePos is a function that detects features in a movie');
                        disp('To work, it needs to receive a structure containing 2 detection parameter:');
                        disp('delta which is the radius around which detection is performed? usually 6 pixels');
                        disp('chi2 which characterize how certain you want to be about the fact that there is a molecule there');
                        disp('Usually between 20 and 200');

                    end
            end
            
            obj.candidatePos = candidate;
            
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
        
        function superResConsolidate(obj,roiSize,frames)
            
            assert(obj.checkStatus('calibrated'),'Data should be calibrated to consolidate');
            assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');    
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before consolidation');
            
            [file2Analyze] = getFileInPath(obj, obj.raw.movInfo.Path, '.mat');
            
            if any(contains({file2Analyze.name},'particle')==true)
                quest = 'Some consolidated positions were found in the raw folder, do you want to load them or run again ?';
                title = 'Question to User';
                btn1  = 'Load';
                btn2 = 'run again';
                defbtn = 'Load';
                answer = questdlg(quest,title,btn1,btn2,defbtn);
                
                switch answer
                case 'Load'
                    particle = load([file2Analyze(1).folder filesep 'particle.mat']);
                    particle = particle.particle;
                    run = false;
                case 'run again'
                    run = true;
                end
            else
                run = true;
            end    
            
            if run
            
                switch nargin

                    case 1

                        frames = 1: obj.calibrated.nFrames;
                        roiSize = 6;
                        disp('Running consolidation on every frame with roi of 6 pixel');

                    case 2

                        frames = 1: obj.calibrated.nFrames;
                        disp('Running consolidation on every frame')

                    case 3

                        [frames] = obj.checkFrame(frames);
                        assert(min(size(roiSize))==1,'RoiSize is expected to be a single number')
                        assert (isnumeric(roiSize),'RoiSize is expected to be a single number');

                    otherwise

                        error('Something wrong with input'); 

                end

                nFrames = length(frames);
                %allocate for storage
                particleList = cell(1,obj.raw.maxFrame(1));
                nParticles = zeros(1,obj.raw.maxFrame(1));
                idx2TP = zeros(1,obj.raw.maxFrame(1));
                h = waitbar(0,'Consolidating candidate ...');
                for i = 1 : 1:nFrames
                    disp(['Consolidating frame ' num2str(i) ' / ' num2str(nFrames)]);
                    idx = frames(i);
                    [data] = obj.getFrame(idx);
                    [candidate] = obj.getCandidatePos(idx);

                    if isempty(candidate)

                        warning('Frame %d did not contain any candidate',idx);
                        particleList{idx} = nan(5);
                        nParticles(idx) = 0;


                    else

                        [finalCandidate] = obj.consolidatePos(data, candidate, roiSize);
                        particleList{idx} = finalCandidate;
                        nParticles(idx) = length(finalCandidate);
                        if ~isempty(finalCandidate)
                            idx2TP(idx) = idx;
                        end

                    end
                    waitbar(i/nFrames,h,['Consolidating candidate... ' num2str(i) '/' num2str(nFrames) ' done']);
                end
                close(h);
                particle.List   = particleList;
                particle.nParticles = nParticles;
                particle.tPoint = nFrames;
                particle.idx2TP = nonzeros(idx2TP);
                particle.roiSize = roiSize;
                particle.Traces = [];
                particle.nTraces = [];

                fileName = sprintf('%s%sparticle.mat',obj.raw.movInfo.Path,'\');
                save(fileName,'particle');
            end
            
            obj.particles = particle;
        end
        
        function [particle] = getParticles(obj,frames)
            
            [idx] = obj.checkFrame(frames);
            particle = obj.particles.List{idx};
            
            if isempty(particle)
                
                warning('There was no particle found in this frames, check than you ran superResConsolidate method on this frame beforehand');
            
            end
        end
        
        function [traces,counter] = trackInZ(obj)
            assert(obj.checkStatus('calibrated'),'Data should be calibrated to do ZCalibration');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before zCalibration');
            assert(~isempty(obj.particles), 'No particles found, please run superResConsolidate method before doing ZCalibration');
            
            %We copy the List as boolean to keep track of where there are
            %still particles left
            [listBool] = obj.copyList(obj.particles.List,1);
            %We copy as NaN for storage of the traces;
            [traces]   = obj.copyList(obj.particles.List,NaN);
            
            [idx] = obj.pickParticle(listBool);
            counter = 1;
            errCount =1;
            while (idx)
                
            if errCount>100
                warning('While loop ran for unexpectedly longer time');
                break;
                
            end
            
            [listIdx] = obj.connectParticles(obj.particles.List,listBool,idx,'ZCal');
           
            if length(listIdx) < 5 
                
                [listBool] = obj.removeParticles(listBool,listIdx);
                
            else
                
            [traces]  = obj.storeTraces(traces,listIdx,counter);
            
            counter = counter +1;
            
            [listBool] = obj.removeParticles(listBool,listIdx);

            end
            
            [idx] = obj.pickParticle(listBool);
            errCount = errCount +1;
            end
            counter = counter -1;
            obj.particles.Traces = traces;
            obj.particles.nTraces = counter;
      
        end
        
        function [zData] = zCalibrate(obj)
            assert(obj.checkStatus('calibrated'),'Data should be calibrated to do ZCalibration');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before zCalibration');
            assert(~isempty(obj.particles), 'No particles found, please run superResConsolidate method before doing ZCalibration');
            
            if ~isempty(obj.particles.Traces)
                quest = 'Some tracked traces were found in the object, do you want to keep them or run again ?';
                title = 'Question to User';
                btn1  = 'Keep';
                btn2 = 'Run again';
                defbtn = 'Keep';
                answer = questdlg(quest,title,btn1,btn2,defbtn);
                
                switch answer
                    case 'Keep'
                    case 'Run again'
                       
                        [traces, nPart] = obj.trackInZ;
                    otherwise
                        error('Unknown answer to question dialog ');
                        
                end
                
            else
                
                %We first track particle in Z
                [traces, nPart] = obj.trackInZ;
                
            end
            %Transform the traces into calibration data (extracting data
            %per plane
            zCalData = obj.getZCalData(obj.particles.Traces,obj.particles.nTraces);
            
            %synchronize the data
            zSyncCalData = obj.syncZCalData(zCalData);
            
            zCal =  obj.calZCalibration(zSyncCalData);
            
            obj.zCalibration.cal = zCal;
            obj.zCalibration.calData = zCalData;
            obj.zCalibration.syncEllip = zSyncCalData;
            
            zData = obj.zCalibration;
        end
        
        function [traces] = get3DTraces(obj)
            list = obj.particles.List;
            tracesIdx = obj.particles.Traces;
            pxSize = obj.info.pxSize;
            elliptRange = [obj.zCalibration.syncEllip{1,3}(1) obj.zCalibration.syncEllip{1,3}(2)];
            elliptRange = elliptRange(1):0.01:elliptRange(2);
            wRange1 = length(elliptRange(elliptRange<=1));
            wRange2 = length(elliptRange(elliptRange>=1));
            
            weight1 = linspace(1,5,wRange1);
            weight2 = linspace(5,1,wRange2);
            
            finalWeight = [weight1 weight2];
            
            if isempty(obj.zCalibration.cal)
                
                warning('no z calibration detected, only show 2D plot');
                
            end
            traces = zeros(length(list),6,obj.particles.nTraces);
            for i = 1 : length(list)
                if ~isempty(list{i})
                    for j = 1 : length(list{i})
                        
                        currentPart = list{i}{j};
                        
                        %find indices to data in correct ellipticity range
                        id = and(currentPart(:,3)>=elliptRange(1),...
                                       currentPart(:,3)<=elliptRange(end));
                        if all(id==0)
                            
                        else
                                  
                        ellip2Keep = currentPart(id,3);
                        idx = ellip2Keep;
                        for k = 1 :length(ellip2Keep)
                            
                            [~,idx(k)] = min(abs(elliptRange-ellip2Keep(k)));
                            
                            
                        end
                           
                        weight = finalWeight(idx);
                        xAvg = sum(diag(list{i}{j}(id,2)* weight))/sum(weight) * pxSize;
                        yAvg = sum(diag(list{i}{j}(id,1)* weight))/sum(weight) * pxSize;
                       
                        x = list{i}{j}(3,2)* pxSize;
                        y = list{i}{j}(3,1)* pxSize;
                        
                        if ~isempty(obj.zCalibration.cal)
                            
                            [z,zAvg] = obj.getZPosition(list{i}{j},elliptRange,finalWeight);
                            
                        else
                            
                            zAvg = 0;
                            z = 0;
                            
                        end
                        
                        if ~isnan(tracesIdx{i}{j})
                        traces(i,:,tracesIdx{i}{j}) = [x y z xAvg yAvg zAvg];
                        else
                        end
                        
                        end
                    end                    
                end
            end
            obj.particles.traces3D = traces;
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
                plot(colPos(planeIdx==i),rowPos(planeIdx==i),'g+')
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
        
        function showParticles(obj,idx)
            assert(length(idx)==1, 'Only one frame can be displayed at once');
            [idx] = obj.checkFrame(idx);
            % Show Candidate
            obj.showCandidate(idx);
            
            if isempty(obj.particles)
                
                warning('You did not consolidate the candidate yet, please use the consolidate method before showing the particles');
                
            else
                
                if isempty(obj.particles.List(idx))
                    
                    warning('The candidates of the requested frame were not consolidated yet, only showing the candidate');
                    
                else
                    
                    roiSize = obj.particles.roiSize;
                    nParticles = obj.particles.nParticles(idx);
                    h = gcf;
                    nPlanes = obj.calibrated.nPlanes;
                    colors = rand(nParticles,3);
                    
                    for i = 1 : nPlanes
                        subplot(2,nPlanes/2,i)
                        hold on
                        for j = 1 : nParticles
                            currPart = obj.particles.List{idx}{j};
                            if(~isempty(currPart(currPart(:,end) == i)))
                                part2Plot = currPart(currPart(:,end) == i,:);
                                plot(part2Plot(2),part2Plot(1),'o',...
                                    'LineWidth',2, 'MarkerSize',10, 'MarkerEdgeColor',colors(j,:));
                            end
                        end
                        hold off
                    end
                    
                    [frame] = getFrame(obj,idx);
                    assert(isstruct(frame),'Error unknown data format, data should be a struct');
                    for i = 1:nParticles
                        
                        currPart = obj.particles.List{idx}{i};
                        %Remove rows containing NaNs
                        idx2NaN = isnan(currPart(:,1));
                        currPart(idx2NaN,:) = [];
                        planes = currPart(:,end);
                        figure(20+i)
                        hold on
                        for j = 1 : length(planes)
                            jdx = planes(j);
                            currFrame = frame.(sprintf('plane%d',jdx));
                            ROI = EmitterSim.getROI(currPart(j,2), currPart(j,1),...
                         roiSize, size(currFrame,2), size(currFrame,1));
                            subplot(1,length(planes),j)
                            imagesc(currFrame(ROI(3):ROI(4),ROI(1):ROI(2)));
                            title({['Particle ' num2str(i)],[ ' Plane ' num2str(jdx)]});
                            axis image
                            colormap('jet')

                        end
                        hold off
                    end
                end
            end              
        end
        
        function showParticlesTracked(obj,ips)
            assert(~isempty(obj.particles.Traces),'You need to get the traces before displaying them');
            
            nPlanes = obj.calibrated.nPlanes;
            colors = rand(obj.particles.nTraces,3);
            
            for i = 1 : length(obj.particles.Traces)
            nParticles = length(obj.particles.Traces{i});
            obj.showCandidate(i);
            
            h = gcf;
            
                for j = 1 : nPlanes
                    subplot(2,nPlanes/2,j)
                    hold on
                    for k = 1 : nParticles
                        currPart = obj.particles.List{i}{k};
                        labelColor = obj.particles.Traces{i}{k};
                        if ~isnan(labelColor)
                            if(~isempty(currPart(currPart(:,end) == j)))
                                part2Plot = currPart(currPart(:,end) == j,:);
                                plot(part2Plot(2),part2Plot(1),'o',...
                                    'LineWidth',2, 'MarkerSize',10, 'MarkerEdgeColor',colors(labelColor,:));
                                drawnow
                            end
                        end
                    end
                    hold off
                    pause(1/ips);
                end
            end
        end
        
        function showTraces(obj)
            assert(~isempty(obj.particles.Traces),'You need to get the traces before displaying them');
            assert(isfield(obj.particles,'traces3D'),'You need to extract the 3D traces before display');
            traces = obj.particles.traces3D;
            npart = size(traces,3);
            %plot XYZ for every particles
            
            figure()
            
            for i = 1:npart
                data = traces(:,:,i);
                data = data(data(:,1)~=0,:);
                fprintf('std in X from best focus: %0.2f \n',nanmedian(nanstd(data(:,1))));
                fprintf('std in y from best focus: %0.2f \n',nanmedian(nanstd(data(:,2))));
                fprintf('std in Z from best focus: %0.2f \n',nanmedian(nanstd(data(:,3))));
                fprintf('std in X from mean planes: %0.2f \n',nanmedian(nanstd(data(:,4))));
                fprintf('std in y from mean planes: %0.2f \n',nanmedian(nanstd(data(:,5))));
                fprintf('std in Z from mean planes: %0.2f \n',nanmedian(nanstd(data(:,6))));
                
                
                subplot(2,1,1)
                hold on
                plot(data(:,1) - data(1,1));
                plot(data(:,2) - data(1,2));
                plot(data(:,3) - data(1,3));
                title('X Y Z position for different particle in best focus')
                xlabel('Frame')
                ylabel('Position(nm)')
                hold off
                
                subplot(2,1,2)
                hold on
                plot(data(:,4) - data(1,4));
                plot(data(:,5) - data(1,5));
                plot(data(:,6) - data(1,6));
                xlabel('Frame')
                ylabel('Position(nm)')
                title('X Y Z position for different particle mean')
                hold off
            end
           
            
            %plot Euclidian distance
             %Calc euclidian distance
                 
            figure()
            hold on
            for i = 1:npart
                data = traces(:,:,i);
                data = data(data(:,1)~=0,:);
                
                eucl = sqrt((data(:,1)-data(1,1)).^2 + (data(:,2)-data(1,2)).^2 +...
                    (data(:,3)-data(1,3)).^2 );
                medEucl = sqrt((data(:,4)-data(1,4)).^2 + (data(:,5)-data(1,5)).^2 +...
                    (data(:,6)-data(1,6)).^2 );
                
                eucl2D = sqrt((data(:,1)-data(1,1)).^2 + (data(:,2)-data(1,2)).^2);
                medEucl2D = sqrt((data(:,4)-data(1,4)).^2 + (data(:,5)-data(1,5)).^2);
                
                fprintf('std in 2D from best focus: %0.2f \n',nanmedian(nanstd(eucl2D)));
                fprintf('std in 2D from mean of planes: %0.2f \n',nanmedian(nanstd(medEucl2D)));
                fprintf('std in 3D from best focus: %0.2f \n',nanmedian(nanstd(eucl)));
                fprintf('std in 3D from mean of planes: %0.2f \n',nanmedian(nanstd(medEucl)));
                
                subplot(2,2,1)
                hold on
                plot(eucl);
                title({'3D euclidian distance'; 'From best focus'})
                hold off
                
                subplot(2,2,2)
                hold on
                plot(medEucl);
                title({'3D euclidian distance'; 'From median'})
                hold off
                
                subplot(2,2,3)
                hold on
                plot(eucl2D);
                title({'2D euclidian distance'; 'From best focus'})
                hold off
                
                subplot(2,2,4)
                hold on
                plot(medEucl2D);
                title({'2D euclidian distance'; 'From median'})
                hold off
                
                %ylim([-400 400]);
                

            end
            hold off
            
            
        end
        
        function showZCalibration(obj)
            assert(~isempty(obj.zCalibration),'No zCalibration found, please run z calibration before display');
            relZ = obj.calibrated.oRelZPos*1000;%in nm
            figure()
            for i = 1 : length(obj.zCalibration.syncEllip)
                
                z =  obj.zCalibration.syncEllip{i}(:,1)+relZ(i);
                ellip = obj.zCalibration.syncEllip{i}(:,2);
                
                ellip1 = ellip(ellip>=1);
                z1 = z(ellip>=1);
                
                ellip2 = 1./ellip(ellip<=1);
                z2 = z(ellip<=1);
               yAx = 1:0.1:2;
               subplot(2,1,1)
               hold on
               scatter(z1,ellip1)
               plot(ones(length(yAx),1)*relZ(i),yAx)
               title('Elliptiticy Elongated in Y')
               xlabel('zPos (nm)')
               ylabel('Ellipticity (sigY/sigX)')
               ylim([1 2])
               hold off
               
               subplot(2,1,2)
               hold on
               scatter(z2,ellip2)
               plot(ones(length(yAx))*relZ(i),yAx)
               title('Elliptiticy Elongated in X')
               xlabel('zPos (nm)')
               ylabel('Ellipticity (sigX/sigY)')
               ylim([1 2])
               hold off
                
            end
            
            figure()
            scatter(obj.zCalibration.syncEllip{1,2}(:,1),obj.zCalibration.syncEllip{1,2}(:,2));
            xlabel('ZPosition')
            ylabel('Ellipticity')
            title('Ellipticity-Z curve for all the planes superimposed')
            
            figure()
            
            for i = 1 : length(obj.zCalibration.syncEllip)
                
                dataCurrentPlane = obj.zCalibration.syncEllip{i};
                
                [binnedData] = obj.zCalBinning(dataCurrentPlane,length(dataCurrentPlane)/5);
                zVec = min(dataCurrentPlane(:,1)):max(dataCurrentPlane(:,1));
                %retrieving fit to display
                p = obj.zCalibration.cal{i};
                fit = polyval(p,zVec);
                %shifting according to the plane
                zVec = zVec + relZ(i) ;
                dataCurrentPlane(:,1) = dataCurrentPlane(:,1)+ relZ(i); 
                binnedData(:,1) = binnedData(:,1) +relZ(i);
                
                subplot(1,2,1)
                hold on
                scatter(binnedData(:,1),binnedData(:,2))
                plot(zVec,fit)
                title('Binned data fitted with polynomial')

                subplot(1,2,2)
                hold on
                scatter(dataCurrentPlane(:,1),dataCurrentPlane(:,2))
                plot(zVec,fit)
                title('Full data fitted with polynomial')
                
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
            candMet = zeros(length(frameCandidate),6);
            nP = numel(fields(data));
            for j = 1 : nP
                
                planeCandidate = frameCandidate(frameCandidate(:,3)==j,1:2);
                planeData = data.(sprintf('plane%d',j));

                cM = zeros(size(planeCandidate,1),5);
                if(~isempty(planeCandidate))
                    for k = 1 : size(planeCandidate,1)
                        %Get the ROI    
                        [roi_lims] = EmitterSim.getROI(planeCandidate(k,2), planeCandidate(k,1),...
                            delta, size(planeData,2), size(planeData,1));
                        ROI = planeData(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2));
                        %Phasor fitting to get x,y,e
                        [row,col,e] = Localization.phasor(ROI);
                        [LRT,~] = Localization.likelihoodRatioTest(ROI,sig,[row col]);
                        rowPos = planeCandidate(k,1) + row;
                        colPos = planeCandidate(k,2) + col;
                        
                        [grad,~] = imgradient(ROI);
                        Grad = max(max(grad,[],2),[],1);

                        cM(k,:) = [rowPos colPos e LRT Grad];

                    end
                    %candidate Metric: phasor localization Res, Likelihood Res,
                    %indexCandidate in plane, Plane.
                    %candidateMetric: [row col e LRT Grad J]

                    candMet(currentk:currentk+k-1,:) = [cM j*ones(k,1)];
                    currentk = currentk+k;
                end
                
            end
            
            [corrEllip, focusMetric] = Localization.calcFocusMetric(candMet(:,3),candMet(:,4));
            %Kill "bad" PSF from the focus metric so they are not taken
            %later
            %reformating to keep the same format as how the data is saved
            %later
            candMet = [candMet candMet(:,end)];
            candMet(:,6) = candMet(:,5);
            candMet(:,5) = focusMetric;
            candMet(:,4) = [];
            focusMetric((1-corrEllip)>0.3) = NaN;
            
            counter = 1;
            nPart = 0;
            maxIt = size(candMet,1);
            finalCandidate = cell(max(size(find(~isnan(focusMetric)))),1);
            
            while and(~isempty(focusMetric), ~isnan(nanmax(focusMetric)))
                
                if counter> maxIt
                    
                    error('While loop ran for an unexpectedly long time, something might be wrong');
               
                end
              
                %Find candidate in best focus
                [~,idx] = max(focusMetric);
                currentPlane = candMet(idx,end);
              
                %Check which planes are to be check
                planes2Check = currentPlane-2:currentPlane-1;
                planes2Check = planes2Check(planes2Check>0);

                planes2Check = [planes2Check currentPlane+1:currentPlane+2];
                planes2Check = planes2Check(planes2Check<nP+1);
                currentCand = candMet(idx,:);
                direction = 1;
                
                particle = nan(5,6);
                particle(3,:) = currentCand;
                nCheck = length(planes2Check);
                for i = 1:nCheck
                    
                    cand = candMet(candMet(:,end) == planes2Check(i),:);
                    if(planes2Check(i) > currentPlane)
                        direction = -1;
                    end
                    
                    [isPart] = obj.isPartner(currentCand,cand,direction,'plane');
                    if ~all(isPart ==0)
                        id = cand(isPart,end)-currentCand(end);
                        particle(round(length(currentCand)/2)+id,:) = cand(isPart,:);
                    end
                    
                end
                
                 
                
                %Check if the resulting configuration of the plane make
                %sense
                planeConfig = particle(:,end);
                [checkRes] = obj.checkPlaneConfig(planeConfig);
                %Store
                if checkRes
                    
                    nPart = nPart +1;
                    finalCandidate{nPart,1} = particle;
                    
                    %We remove the particle(s) from the list
                    focusMetric(ismember(candMet(:,1), particle(:,1))) = [];
                    candMet(ismember(candMet(:,1), particle(:,1)),:) = [];
                    
                    
                else
                    %Otherwise we remove it from the best focus search list
                    %by putting focus metric to NaN
                    focusMetric(idx) = NaN;
                    
                end

                counter = counter+1;
                
            end
            %we delete empty cells from the array
            idx2Empty = cellfun(@isempty,finalCandidate);
            finalCandidate(idx2Empty(:,1),:) = [];
            
end
        
        function [isPart] = isPartner(obj, current, next, direction, check)
            
            %This function is designed to have PSFE plate ON
            assert(abs(direction) == 1, 'direction is supposed to be either 1 (up) or -1 (down)');
            assert(size(current,2) == size(next,2), 'Dimension mismatch between the tail and partner to track');
           
            switch check
                case 'plane' %TL, consolidation between Plane
                    Thresh = 2*sqrt(2); %in Px As the superResCal is not 
                    %performed yet we allow 2px in both direction
                    [checkRes1] = obj.checkEuDist(current(:,1:2),...
                                next(:,1:2),Thresh);

                    % Test ellipticity
                    [checkRes2] = obj.checkEllipticity(current(:,3),...
                        next(:,3),direction);
                    
                    % Test focus Metric    
                     maxExpFM = current(4)+0.1*current(4);
                     checkRes3 = next(:,4) < maxExpFM;  
                           
                   %isPart will only be true for particle that passes the 3 tests       
                   isPart = checkRes1.*checkRes2.*checkRes3;
                   
                   if(length(find(isPart))>1)

                        warning('Could not choose which particle was the partner of the requested particle, killed them both');
                        isPart(isPart==1) = 0;

                   end

                case 'ZCal' %ZStack, consolidation between frame
                    %The calculation here is ran in parallel, we check if
                    %the current particle is a partner of one of the
                    %particles in the next frame. Some indexing or step
                    %might therefore seems unnecessary but allow to take
                    %any number of particles
                    
                    isEdgePlane = or(~isempty(find(current(:,end)==1,1)),~isempty(find(current(:,end)==8,1)));
                    
                    if isEdgePlane
                        
                        nConsistentPlanes = 2;
                        
                    else
                        
                        nConsistentPlanes = 2;
                         
                    end
                    
                    
                    isPart = zeros(size(next,3),1);
                    %check focus is not more than one plane away 
                    roughcheck1 = squeeze(abs(current(3,end)-next(3,end,:)<=1));
                    
                    if all(roughcheck1 ==0)
                        disp('Something is wrong your focus changed more than one plane between 2 frames');
                    else
                        %Check that at least 2 planes are in common
                        commonPlanes = obj.findCommonPlanes(current(:,end),squeeze(next(:,end,:)));
                        roughcheck2  = sum(commonPlanes,1)>=nConsistentPlanes;
                        if all(roughcheck2 ==0)
                            disp('Less than 2 planes in common, breaking out');
                        else 
                            
                            for i = 1 : size(next,3)
                                % Test Euclidian distance
                                Thresh = 1; %in px
                                [checkRes1] = obj.checkEuDist(current(commonPlanes(:,1,i),1:2),...
                                    squeeze(next(commonPlanes(:,2,i),1:2,i)),Thresh);

                                % Test ellipticity
                                eWeight = [1 2 3 2 1];
                                thresh = 5;
                                [checkRes2] = obj.checkEllipticity(current(commonPlanes(:,1,i),3),...
                                    squeeze(next(commonPlanes(:,2,i),3,i)),direction,thresh,eWeight(commonPlanes(:,1,i)));
                                 %To be a particle, we want the position to be
                                %consistent in at least 2 planes and
                                %ellipticity to pass the test.
                                isPart(i) = and(length(find(checkRes1))>=nConsistentPlanes, checkRes2);
                                
                            end
                                              
                        end
                    end
                              
                                
                    
                case 'Track' %TL, consolidation between frame
                otherwise 
                    error('Unknown type of experiment');
            end
            isPart = logical(isPart);
         
                
        end
     
        function [commonPlanes] = findCommonPlanes(~,planeInCurrent,planeInNext)
            
            commonPlanes = zeros(size(planeInNext,1),2,size(planeInNext,2));
            
            for i = 1 : size(planeInNext,2)
                
                commonPlanes(:,1,i) = ismember(planeInCurrent,planeInNext(:,i));
                commonPlanes(:,2,i) = ismember(planeInNext(:,i),planeInCurrent);
                
            end
            commonPlanes = logical(squeeze(commonPlanes));
        end
        
        function listIdx = connectParticles(obj,List,listBool,idx,criteria)

            isPart = true;
            counter = 1;

            listIdx = zeros(length(List)-idx(1),2);
            listIdx(1,:) = idx;
            currentIdx = idx;
            while isPart
                if currentIdx >= length(List)
                           break;
                end
                part2Track = List{currentIdx(1)}{currentIdx(2)};
                [checkRes] = obj.checkListBool(listBool,currentIdx(1)+1);
                nPartInFrame = length(checkRes);

                if ~all(checkRes==0)
                   %We use reshape to input the different particles of the next
                   %frame at once by storing them in the 3rd dimension
                   shape = [size(part2Track,1) size(part2Track,2) nPartInFrame]; 
                   nextParts(:,:,1:nPartInFrame) = reshape(cell2mat(List{currentIdx(1)+1}'),shape);
                   nextParts = nextParts(:,:,logical(checkRes));
                   [isPart] = obj.isPartner(part2Track,nextParts,1,criteria);

                   counter = counter+1;
                   if(length(find(isPart==1))>1)

                       warning('Could not choose between 2 close particles, killing them both')
                       isPart = false;

                   elseif (~all(isPart==0))
                       %Update newIdx
                       currentIdx(1) = currentIdx(1)+1;%current become the connected from next frame
                       idx2AvailableParticles = find(checkRes);
                       currentIdx(2) = idx2AvailableParticles(isPart);
                       nextParts = [];
                       listIdx(counter,:) = currentIdx;

                       isPart = true;
                   else

                       isPart = false;

                   end

                   if counter == length(obj.particles.List)-1

                       isPart = false;

                   end
                   
                else
                    isPart = false;
                end
            listIdx(listIdx(:,1) == 0,:) = [];

            end
        end
        
        function [listCopy] = copyList(~, List, filling)
            
            assert(iscell(List),'The List should be a cell array');
            listCopy = cell(List);
            
            for i = 1:length(List)
                
                if isempty(List{i})
                    
                    listCopy{i} = [];
                    
                else
                    
                    for j = 1:length(List{i})
                        
                        listCopy{i}{j} = filling;
                        
                    end
                    
                end
                
            end
            
        end
        
        function [idx] = pickParticle(obj,listBool)
            %This function will pick a particle and return false if there
            %is none
              
            for i = 1:length(listBool)

                if exist('idx','var')

                    break;

                else

                    if isempty(listBool{i})

                    else

                        for j = 1:length(listBool{i})

                            if listBool{i}{j} == true

                                idx = [i,j];
                                break;

                            end          
                        end
                    end
                end
            end
              
            if ~exist('idx','var')

              idx = false;

            end
            
        end
        
        function [traces] = storeTraces(obj,traces,listIdx,numPart)
            
            for i = 1 : size(listIdx,1)
                
                traces{listIdx(i,1)}{listIdx(i,2)} = numPart;
                
            end
            
        end
        
        function [listBool] = removeParticles(obj,listBool,listIdx)
            
            for i = 1 : size(listIdx,1)
                
                listBool{listIdx(i,1)}{listIdx(i,2)} = false;
                
            end
            
        end
        
        function [zCalData] = getZCalData(obj,traces,nPart)
            
            
            zCalData = cell(obj.calibrated.nPlanes,nPart);
            
               for i = 1:length(traces)

                    if isempty(traces{i})

                    else

                        for j = 1:length(traces{i})
                            if ~isnan(traces{i}{j})
                            planes = obj.particles.List{i}{j}(:,end);
                            planes = planes (~isnan(planes));
                            
                            for k = 1 : length(planes)
                                data2Store = obj.particles.List{i}{j};
                                data2Store = data2Store(data2Store(:,end) == planes(k),:);
                                zCalData{planes(k),traces{i}{j}}(i,:) = [data2Store(3:end) i] ;
                                
                            end
                            else
                            end
                            
                        end
                    end
               end
               %clean up the data
               for i = 1 : size(zCalData,1)
                   for j = 1 : size(zCalData,2)
                       
                       if ~isempty(zCalData{i,j})
                           
                            zCalData{i,j}(zCalData{i,j}(:,1)==0,:) = [];
                            
                       else
                       end
                       
                   end
               end
        end   
        
        function [zSyncCalData] = syncZCalData(obj,zCalData)
            deg = 4;
            minEllipt = 0.77;
            maxEllipt = 1.6;
            zStep = obj.getZStep;
            zSyncCalData = cell(8,2);
%             figure
            for i = 1:size(zCalData,1)
                for j = 1:size(zCalData,2)
                    %Extrac data in acceptable ellipticity range
                    %(arbitrary)
                    if ~isempty(zCalData{i,j})
                        zPos = zCalData{i,j}(:,5)*zStep;
                        ellipt = zCalData{i,j}(:,1);
                        zPos = zPos(and(ellipt<maxEllipt, ellipt > minEllipt));
                        ellipt = ellipt(and(ellipt<maxEllipt, ellipt > minEllipt));

                        %Now we shift in z the value closest to ellipt =1 for
                        %the fit == rough synchronization    
                        [~,idx] = min(abs(ellipt-1));
                        focus1 = zPos(idx);
                        zPos = zPos - focus1;

                        if length(ellipt)<deg
                        else
                            %Do the fit and extract the z position of the focus
                            p = polyfit(zPos,ellipt,deg);
                            zVec = min(zPos):0.001:max(zPos);%1nm step

                            if or(min(zPos)>-0.1, max(zPos) <0.1)
                                zVec = -1:0.001:1; 
                            end

                            fit = polyval(p,zVec);
                            [~,idx] = min(abs(fit-1));
                            focus2 = zVec(idx);                  
                            shift = focus1+focus2;
                            zCalData{i,j}(:,6) = ((zCalData{i,j}(:,5))*zStep -shift)*1000;
                            fullSync = ((zCalData{i,j}(:,5))*zStep -shift)*1000;
                            zSyncCalData{i,1} = [zSyncCalData{i,1}; zCalData{i,j}(:,[6 1])];
                            zSyncCalData{1,2} = [zSyncCalData{1,2}; [fullSync zCalData{i,j}(:,1)]];
                        end
                    end
                end
                
                [~,ind] = sort(zSyncCalData{i,1}(:,1));
                zSyncCalData{i,1} = zSyncCalData{i,1}(ind,:);
                
            end
            
            [~,ind] = sort(zSyncCalData{1,2}(:,1));
            zSyncCalData{1,2} = zSyncCalData{1,2}(ind,:);
            zSyncCalData{1,3} = [minEllipt, maxEllipt, deg];
            
        end
        
        function [zStep] = getZStep(obj)
            nFrames = obj.raw.maxFrame(1);
            zPos = zeros(nFrames,1);
            
            for i = 1 : obj.raw.maxFrame(1)
                
                zPos(i) = obj.raw.frameInfo(2*i).Pos(3);
                
            end
            
            zStep = mean(diff(zPos));
            
        end
        
        function [zCalibration] = calZCalibration(obj, zSyncCalData)
            zCalibration = cell(length(zSyncCalData),1);
            deg = zSyncCalData{1,3}(3);
            
            for i = 1: length(zSyncCalData)
                dataCurrentPlane = zSyncCalData{i};
                [binnedData] = obj.zCalBinning(dataCurrentPlane,length(dataCurrentPlane)/5);
                
                zVec = min(dataCurrentPlane(:,1)):max(dataCurrentPlane(:,1));
     
                p = polyfit(dataCurrentPlane(:,1),dataCurrentPlane(:,2),deg);
                fit = polyval(p,zVec);
%                 
%                 figure
%                 subplot(1,2,1)
%                 hold on
%                 scatter(binnedData(:,1),binnedData(:,2))
%                 plot(zVec,fit)
%                 
%                 subplot(1,2,2)
%                 hold on
%                 scatter(dataCurrentPlane(:,1),dataCurrentPlane(:,2))
%                 plot(zVec,fit)
 
                zCalibration{i} = p;
                
            end
            
        end
        
        function [binnedData] = zCalBinning(obj,zCalData2Bin,nPoint)
            
            nQuant = linspace(0,1,length(zCalData2Bin)/nPoint);
            binnedData = quantile(zCalData2Bin,nQuant);
            
        end
        
        function [zPos,zAvg] = getZPosition(obj,particle,elliptRange,finalWeight)
             assert(~isempty(obj.zCalibration),'No zCalibration found, please run z calibration calculating Z position');
             zCal = obj.zCalibration.cal;
             relZ = obj.calibrated.oRelZPos;
             syncEllip = obj.zCalibration.syncEllip;
             
             zVec = min(syncEllip{particle(3,end),1}) : 1 : max(syncEllip{particle(3,end),1});
             fit = polyval(zCal{particle(3,end),1},zVec);
             %find the index of the value the closest to the particle
             %ellipticity
             [~,idx] = min(abs(fit-particle(3,3)));
            zPos = zVec(idx)+ relZ(particle(3,end))*1000;
             
             
            id = and(particle(:,3)>=elliptRange(1),...
                                       particle(:,3)<=elliptRange(end));
            data2Keep = particle(id,:);
            
           
             
            idx1 = zeros(size(data2Keep,1),1);
            zAvg = zeros(size(data2Keep,1),1);
           
            for k = 1 :size(data2Keep,1)

                [~,idx1(k)] = min(abs(elliptRange-data2Keep(k,3)));
                
                zVectmp = min(syncEllip{data2Keep(k,end),1}) : 1 : max(syncEllip{data2Keep(k,end),1});
                fit = polyval(zCal{data2Keep(k,end),1},zVectmp);
                
                [~,idx] = min(abs(fit-data2Keep(k,3)));
                zAvg(k) = zVectmp(idx) + relZ(data2Keep(k,end))*1000;
                

            end
                     
             weight = finalWeight(idx1);               
             zAvg = sum(diag(zAvg(:)* weight(:)'))/sum(weight);

             
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
        
        function [frames]   = checkFrame(obj, frames) 
            
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
        
        function [checkRes] = checkEuDist(obj,current,next,Thresh)

            EuDist = sqrt((current(:,1) - next(:,1)).^2 +...
                            (current(:,2) - next(:,2)).^2);
            checkRes = EuDist < Thresh;
            
            if isempty(checkRes)
                checkRes = false;
            end
            
        end
        
        function [checkRes] = checkEllipticity(obj, current, next, direction, thresh, eWeight)
            
            
            switch direction

            case 1
                
                ellip = current < next+0.1*next;

            case -1

                ellip = current +0.1 *current > next;

            end
            
            if size(current,1)>1
                
                checkRes = sum(ellip .* eWeight(:))>=thresh;
                
            else
                
                checkRes = ellip;
                
            end
            
            if isempty(checkRes)
                checkRes = false;
            end
                    
        end
        
        function [checkRes] = checkPlaneConfig(obj, planeConfig)
            %Here we will check that the consolidation found based on the
            %best focused particle make sense with what we would expect and
            %also that we have enough planes.
            assert(length(planeConfig) <= obj.calibrated.nPlanes,'There is something wrong with your consolidated index and your candidate plane List');
            %Let us test that we have consolidate the particle in at least
            %3 Planes
            isEdgePlane = or(~isempty(find(planeConfig==1,1)),~isempty(find(planeConfig==8,1)));
            
            if isEdgePlane
                
                testPlanes = length(find(~isnan(planeConfig)==true)) >= 2;
                
            else
                
                 testPlanes = length(find(~isnan(planeConfig)==true)) >= 3;
                
            end
           
            
            if testPlanes
               %We check that there is no "Gap" in the plane configuration
               %as it would not make sense.
               testConsec = diff(planeConfig(~isnan(planeConfig)));
               checkRes = all(testConsec==1);
               
            else
                
                checkRes = false;
                
            end
        end
        
        function [checkRes] = checkListBool(~, listBool, idx)
            
            if( idx>= length(listBool))
                
                checkRes = false;
                
            else
            
            if isempty(listBool{idx})
                
                checkRes = false;
                
            else
                test = zeros(1,length(listBool{idx}));
                
                for i = 1 : length(listBool{idx})
                    
                    test(i) = listBool{idx}{i};
                    
                end
                
                if(all(test == 0))
                    
                    checkRes = false;
                    
                else
                    
                    checkRes = test;
                                      
                end 
                
            end
            
            end
        end
    end
end
