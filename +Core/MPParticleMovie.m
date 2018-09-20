classdef MPParticleMovie < Core.MPMovie
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'protected')
        
        candidatePos
        unCorrLocPos
        corrLocPos
        particles
        
    end
    
    methods
        function obj = MPParticleMovie(raw,cal)
            
            obj  = obj@Core.MPMovie(raw,cal);
            
        end
        
        function set.candidatePos(obj,candidatePos)
            
            assert(iscell(candidatePos), 'Expecting a cell containing candidate positions');
            obj.candidatePos = candidatePos;
            
        end
        
        function findCandidatePos(obj,detectParam, frames)
            %Method to perform localization on each plane for each frame
            %Check if some candidate exists already in the folder (previously saved)
            [run, candidate] = Core.MPParticleMovie.existCandidate(obj.raw.movInfo.Path, '.mat');
            
            if run
                switch nargin
                    case 2
                        
                        frames = 1: obj.calibrated.nFrames;
                        disp('Running detection on every frame');
                        
                    case 3
                        
                        [frames] = obj.checkFrame(frames);
                        
                    otherwise
                        
                        error('too many inputs');
                        
                end
                %Localization occurs here
                assert(~isempty(obj.info), 'Missing information about setup to be able to find candidates, please use giveInfo method first or load previous data');
                assert(nargin>1,'not enough input argument or accept loading of previous data (if possible)');
                [candidate] = obj.detectCandidate(detectParam,frames);
                
            elseif ~isempty(candidate)
            else
                %help message
                disp('getCandidatePos is a function that detects features in a movie');
                disp('To work, it needs to receive a structure containing 2 detection parameter:');
                disp('delta which is the radius around which detection is performed usually 6 pixels');
                disp('chi2 which characterize how certain you want to be about the fact that there is a molecule there');
                disp('Typically between 20 and 200');
                
            end
            
            %save the data
            fileName = sprintf('%s%scandidatePos.mat',obj.raw.movInfo.Path,'\');
            save(fileName,'candidate');
            
            obj.candidatePos = candidate;
            obj.info.detectParam = detectParam;
        end
        
        function [candidate] = getCandidatePos(obj, frames)
            %Extract the position of the candidate of a given frame
            [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
            candidate = obj.candidatePos{idx};
            
            if isempty(candidate)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                
            end
        end  
        
        function SRLocalizeCandidate(obj,roiSize,frames)
            assert(~isempty(obj.calibrated),'Data should be calibrated to consolidate');
            assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before consolidation');
            
            [run,locPos] = obj.existLocPos(obj.raw.movInfo.Path,'.mat');
            
            if run
                switch nargin

                    case 1
                        roiSize = 6;
                        frames = 1: obj.calibrated.nFrames;
                        disp('Running SRLocalization on every frame with ROI of 6 pixel radius');

                    case 2

                        frames = 1: obj.calibrated.nFrames;
                        disp('Running SRLocalization on every frame');

                    case 3

                        [frames] = obj.checkFrame(frames);

                    otherwise

                        error('too many inputs');

                end
               
                locPos = cell(size(obj.candidatePos));
                h = waitbar(0,'Fitting candidates ...');
                nFrames = length(frames);
                %Localization occurs here
                for i = 1 : 1:nFrames
                    disp(['Fitting candidates: frame ' num2str(i) ' / ' num2str(nFrames)]);
                    idx = frames(i);
                    %#1 Extract Candidate Position for specific frame
                    [data] = obj.getFrame(idx);
                    [frameCandidate] = obj.getCandidatePos(idx);
                    
                    if isempty(frameCandidate)
                        
                        warning('Frame %d did not contain any candidate',idx);
                        locPos{i} = [];
                        
                    else
                        
                        locPos{i} = obj.superResLocFit(data,frameCandidate,roiSize);
                        
                    end
                    waitbar(i/nFrames,h,['Fitting candidates: frame ' num2str(i) '/' num2str(nFrames) ' done']);
                end
                close(h)
            else
            end
                %save the data
            fileName = sprintf('%s%sSRLocPos.mat',obj.raw.movInfo.Path,'\');
            save(fileName,'locPos');
            
            %store in the object
            obj.unCorrLocPos = locPos;
            obj.corrLocPos   = locPos;
        end
        
        function [locPos] = getLocPos(obj,frames)
             %Extract the position of the candidate of a given frame
            [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
            locPos = obj.unCorrLocPos{idx};
            
            if isempty(locPos)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                
            end
        end
        
        function consolidatePlanes(obj,roiSize,frames)
            %Consolidation refers to connect molecules that were localized
            %at similar position in different plane on a single frame.
            assert(~isempty(obj.calibrated),'Data should be calibrated to consolidate');
            assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before consolidation');
            assert(~isempty(obj.unCorrLocPos),'Localization needs to be performed before consolidation');
           
            %Check if some particles were saved already.
            [run, particle] = Core.MPParticleMovie.existParticles(obj.raw.movInfo.Path, '.mat');
            
            if run
                %Check the number of function input
                switch nargin
                    case 1
                        
                        frames = 1: obj.calibrated.nFrames;
                        roiSize = 6;
                        disp('Running consolidation on every frame with roi of 6 pixel');
                        
                    case 2
                        
                        frames = 1: obj.calibrated.nFrames;
                        disp('Running consolidation on every frame')
                        
                    case 3
                        
                        [frames] = Movie.checkFrame(frames,obj.raw.maxFrame(1));
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
                
                %Consolidation occurs here
                for i = 1 : 1:nFrames
                    disp(['Consolidating frame ' num2str(i) ' / ' num2str(nFrames)]);
                    idx = frames(i);
                    %#1 Extract localized Position for specific frame
                    [fCandMet] = obj.getLocPos(idx);
                    
                    if isempty(fCandMet)
                        
                        warning('Frame %d did not contain any localized positions',idx);
                        particleList{idx} = [];
                        %particleList{idx}{1} = nan(5);
                        nParticles(idx) = 0;
                        
                    else
                        %#2 Consolidate the position of the given frame
                        %across plane
                          %Calculate a focus metric (FM) combining ellipticity and GLRT FM.
                            [corrEllip, focusMetric] = Localization.calcFocusMetric(fCandMet.ellip,fCandMet.fMetric);

                            %reformating to keep the same format as how the data is saved
                            %later
                            fCandMet.fMetric = focusMetric;
                            
                            %focusMetric((1-corrEllip)>0.3) = NaN;
                            
                            %Plane Consolidation occur here
                            [part] = obj.planeConsolidation(fCandMet,focusMetric);

                            %we delete empty cells from the array
                            idx2Empty = cellfun(@isempty,part);
                            part(idx2Empty(:,1),:) = [];
                   
                            particleList{idx} = part;
                            nParticles(idx) = length(part);
                            
                            if ~isempty(part)
                                idx2TP(idx) = idx;
                            end
                        
                    end
                    waitbar(i/nFrames,h,['Consolidating candidate... ' num2str(i) '/' num2str(nFrames) ' done']);
                end
                close(h);
                
                %#3 Storing List
                particle.List       = particleList;
                particle.nParticles = nParticles;
                particle.tPoint     = nFrames;
                particle.idx2TP     = nonzeros(idx2TP);
                particle.roiSize    = roiSize;
                particle.Traces     = [];
                particle.nTraces    = [];
                
                fileName = sprintf('%s%sparticle.mat',obj.raw.movInfo.Path,'\');
                save(fileName,'particle');
            end
            %#4 Storing particles in the object
            obj.particles = particle;
        end
        
        function [particle] = getParticles(obj,frames)
            %GetParticles
            [idx] = obj.checkFrame(frames);
            particle = obj.particles.List{idx};
            
            if isempty(particle)
                
                warning('There was no particle found in this frames, check than you ran superResConsolidate method on this frame beforehand');
                
            end
        end
        
        function showCandidate(obj,idx)
            %Display Candidate
            assert(length(idx)==1, 'Only one frame can be displayed at once');
            [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
            assert(~isempty(obj.candidatePos{idx}),'There is no candidate found in that frame, check that you ran the detection for that frame');
            
            [frame] = getFrame(obj,idx);
            assert(isstruct(frame),'Error unknown data format, data should be a struct, wrong output from getFrame');
            
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
            %display particles (after consolidation), On top of the
            %localization, consolidated particles are circled.
            assert(length(idx)==1, 'Only one frame can be displayed at once');
            [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
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
                    %Display circled
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
                    %Here we display a zoom onto the particle visible on
                    %the specific frame onto the consolidated planes
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
        
        function [candidateList] = planeConsolidation(obj,candMet,focusMetric)
            %Loop through all candidate of a given frame and match them
            %between frame until none can be match or all are matched.
            nPlanes = obj.calibrated.nPlanes;
            counter = 1;
            nPart = 0;
            maxIt = size(candMet,1);
            candidateList = cell(max(size(find(~isnan(focusMetric)))),1);
            %continue until the list is not empty
            while and(~isempty(focusMetric), ~isnan(nanmax(focusMetric)))
                
                if counter> maxIt
                    
                    error('While loop ran for an unexpectedly long time, something might be wrong');
                    
                end
                
                %Find candidate in best focus
                [~,idx] = max(focusMetric);
                currentPlane = candMet.plane(idx);
                
                %Check which planes are to be checked (currently 2 planes
                %above and 2 planes below the given plane
                planes2Check = currentPlane-2:currentPlane-1;
                planes2Check = planes2Check(planes2Check>0);
                planes2Check = [planes2Check currentPlane+1:currentPlane+2];
                planes2Check = planes2Check(planes2Check<nPlanes+1);
                currentCand = candMet(idx,:);
                direction = -1;%Start by checking above
                
                particle = array2table(nan(5,size(currentCand,2)));
                particle.Properties.VariableNames = currentCand.Properties.VariableNames;
                
                particle(3,:) = currentCand;
                nCheck = length(planes2Check);
                camConfig = obj.calibrated.camConfig;
                for i = 1:nCheck
                    
                    cand = candMet(candMet.plane == planes2Check(i),:);
                    if(planes2Check(i) > currentPlane)
                        direction = +1;%check below (Plane 1 is the uppest plane 8 is lowest)
                    end
                    
                    [isPart] = Core.MPParticleMovie.isPartPlane(currentCand,cand,direction);
                    if ~all(isPart ==0)
                        id = cand.plane(isPart)-currentCand.plane;
                        particle(3+id,:) = cand(isPart,:);
                    end
                    
                end
                
                %Check if the resulting configuration of the plane make
                %sense e.g. no hole in the configuration
                
                planeConfig = particle.plane;
                
                [checkRes] = Core.MPParticleMovie.checkPlaneConfig(planeConfig,nPlanes,camConfig);
                %Store
                if checkRes
                    
                    nPart = nPart +1;
                    candidateList{nPart,1} = particle;
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
        end
                

    end
     methods (Static)
        %method linked to candidate
        function [run,candidate] = existCandidate(Path,ext)
            
            [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
            
            %Check if some candidate were already stored
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
                        candidate =[];
                        
                    otherwise
                        error('Unknown answer to user input dialog, most likely due to cancelation')
                end
                
            else
                
                run = true;
                candidate =[];
            end
        end
        
        %method linked to fitting
        function [run,SRLocPos] = existLocPos(Path,ext)
            
            [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
            
            %Check if some candidate were already stored
            if any(contains({file2Analyze.name},'SRLocPos')==true)
                quest = 'Some fitted positions were found in the raw folder, do you want to load them or run again ?';
                title = 'Question to User';
                btn1  = 'Load';
                btn2 = 'run again';
                defbtn = 'Load';
                answer = questdlg(quest,title,btn1,btn2,defbtn);
                
                switch answer
                    case 'Load'
                        
                        SRLocPos = load([file2Analyze(1).folder filesep 'SRLocPos.mat']);
                        name = fieldnames(SRLocPos);
                        SRLocPos = SRLocPos.(name{1});
                        run = false;
                        
                    case 'run again'
                        
                        run = true;
                        SRLocPos =[];
                        
                    otherwise
                        error('Unknown answer to user input dialog, most likely due to cancelation')
                end
                
            else
                
                run = true;
                SRLocPos =[];
            end
        end
        
        %method Linked to particles/planeConsolidation
        function [run, particle] = existParticles(Path, ext)
            
            [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
            %Check if some particles were already saved in the raw folder.
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
                        particle = [];
                        
                    otherwise
                        
                        error('Unknown answer to user input dialog, most likely due to cancelation')
                        
                end
            else
                run = true;
                particle = [];
            end
            
        end
        
        function [isPart]   = isPartPlane(current, next, direction)
            %This function aim at determining whether a candidate from one
            %plane and the another are actually the same candidate on
            %different plane or different candidate. The decision is based
            %on threshold on localization distance, ellipticity and focus
            %metric.
            
            %This function is designed to have PSFE plate ON
            assert(abs(direction) == 1, 'direction is supposed to be either 1 (up) or -1 (down)');
            assert(size(current,2) == size(next,2), 'Dimension mismatch between the tail and partner to track');

            thresh = 4;
            [checkRes1] = Core.MPParticleMovie.checkEuDist([current.row, current.col],...
                [next.row, next.col],thresh);
            
%             % Test ellipticity
            [checkRes2] = Core.MPParticleMovie.checkEllipticity(current.ellip,...
                next.ellip,direction);
            
            % Test focus Metric
            maxExpFM = current.fMetric+0.1*current.fMetric;
            checkRes3 = next.fMetric < maxExpFM;
            
            %isPart will only be true for particle that passes the 3 tests
            isPart = checkRes1.*checkRes2.*checkRes3;
            
            if(length(find(isPart))>1)
                
                warning('Could not choose which particle was the partner of the requested particle, killed them both');
                isPart(isPart==1) = 0;
            end
            
            isPart = logical(isPart);
            
        end
        
        function [checkRes] = checkEuDist(current,next,Thresh)
            %Use to check if the Euclidian distance is within reasonable
            %range
            EuDist = sqrt((current(:,1) - next(:,1)).^2 +...
                (current(:,2) - next(:,2)).^2);
            checkRes = EuDist < Thresh;
            
            if isempty(checkRes)
                checkRes = false;
            end
            
        end
        
        function [checkRes] = checkEllipticity(current, next, direction)
            %Use to check if the Ellipticity make sense with what we
            %expect from the behavior of the PSFEngineering plate
            switch direction
                
                case 1
                    
                    ellip = current < next+0.1*next;
                    
                case -1
                    
                    ellip = current +0.1 *current > next;
                    
            end
            
            checkRes = ellip;

            if isempty(checkRes)
                checkRes = false;
            end
            
          
            
        end
        
        function [checkRes] = checkPlaneConfig(planeConfig,nPlanes,camConfig)
            %Here we will check that the consolidation found based on the
            %best focused particle make sense with what we would expect and
            %also that we have enough planes.
            assert(length(planeConfig) <= nPlanes,'There is something wrong with your consolidated index and your candidate plane List');
            %Let us test that we have consolidate the particle in at least
            %3 Planes
            isEdgePlane = or(~isempty(find(planeConfig==1,1)),~isempty(find(planeConfig==8,1)));
            
            if isEdgePlane
                
                testPlanes = length(find(~isnan(planeConfig)==true)) >= 2;
                
            else
                switch camConfig
                    case 'fullRange'
                        testPlanes = length(find(~isnan(planeConfig)==true)) >= 2;
                    case 'alternated'
                        testPlanes = length(find(~isnan(planeConfig)==true)) >= 3;
                    otherwise
                        error('unknown camera configuration');
                end
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
        
       
     end
     
     methods (Access = private)
         
        %Methods linked to Candidate
        function [candidate] = detectCandidate(obj,detectParam,frames)
            %Do the actual localization
            assert(~isempty(obj.calibrated),'Data should be calibrated to detect candidate');
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
                
                position = table(zeros(500,1),zeros(500,1),zeros(500,1),...
                    zeros(500,1),'VariableNames',{'row', 'col', 'meanFAR','plane'});
                [volIm] = obj.getFrame(frames(i));
                nameFields = fieldnames(volIm);
                
                for j = 1:length(nameFields)
                    %localization occurs here
                    [ pos, meanFAR, ~ ] = Localization.smDetection( double(volIm.(nameFields{j})),...
                        delta, FWHM_pix, chi2 );
                    if ~isempty(pos)
                        startIdx = find(position.row==0,1,'First');
                        pos(:,3) = meanFAR;
                        pos(:,4) = j;
                        position(startIdx:startIdx+size(pos,1)-1,:) = array2table(pos);
                    else
                    end
                end
                
                idx = find(position.row==0,1,'First');
                if isempty(idx)
                    
                    candidate{frames(i)} = position;
                    
                else
                    
                    candidate{frames(i)} = position(1:idx-1,:);
                    
                end
                waitbar(i/nFrames,h,...
                    sprintf('detection of candidates in Frame %d/%d done',i,nFrames));
            end
            
            close(h);
        end
                
        function [candMet] = superResLocFit(obj,data,frameCandidate,roiSize)
            %Candidate metric are determined here (x,y,e,+focusmetric)
            delta = roiSize;
            
            %initialize table
            varNames = {'row','col','z','ellip','magX','magY','meanFAR','fMetric','gFitMet','plane'};
                candMet = table(zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    zeros(size(frameCandidate,1),1),zeros(size(frameCandidate,1),1),...
                    'VariableNames',varNames);
                
            for i = 1:size(frameCandidate,1)
                
                plane = frameCandidate.plane(i);
                planeData = data.(sprintf('plane%d',plane));
                sig = [obj.info.sigma_px obj.info.sigma_px];
                %Get the ROI
                [roi_lims] = EmitterSim.getROI(frameCandidate.col(i), frameCandidate.row(i),...
                    delta, size(planeData,2), size(planeData,1));
                ROI = planeData(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2));
                %Phasor fitting to get x,y,e
                [row,col,e,magX,magY] = Localization.phasor(ROI);
                
                 %LRT focus metric
                [fMetric,~] = Localization.likelihoodRatioTest(ROI,sig,[row col]);
                
                if magX>=magY
                    sig(1) = sig(1) * magX/magY;
                else
                    sig(2) = sig(2) * magY/magX;
                end
                
                %LRT focus metric
                [gFitMet,~] = Localization.likelihoodRatioTest(ROI,sig,[row col]);
                
                rowPos = frameCandidate.row(i) + row;
                colPos = frameCandidate.col(i) + col;

%                 [grad,~] = imgradient(ROI);
%                 Grad = max(max(grad,[],2),[],1);
                
                %storing info
                candMet.row(i) = rowPos;
                candMet.col(i) = colPos;
                candMet.z(i) = 0;
                candMet.ellip(i) = e;
                candMet.magX(i) = magX;
                candMet.magY(i) = magY;
                candMet.meanFAR(i) = frameCandidate.meanFAR(i);
                candMet.fMetric(i) = fMetric;
                candMet.gFitMet(i) = gFitMet;
                candMet.plane(i) = plane;
            end

        end
              
     end
end

