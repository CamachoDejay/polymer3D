classdef MPLocMovie < Core.MPMovie
    %MPLocMovie hold all the method and info necessary for localization,
    %expension of this object will be first centered around Tracking but
    %it will also be easy to extend to STORM/PALM.
    
    properties (SetAccess = 'private')
        
        candidatePos
        particles
        
    end
    
    methods
        
        function obj = MPLocMovie(raw, cal)
            
            obj  = obj@Core.MPMovie(raw,cal);
            
        end
        
        function set.candidatePos(obj,candidatePos)
            
            assert(iscell(candidatePos), 'Expecting a cell containing candidate positions');
            obj.candidatePos = candidatePos;
            
        end
        
        function findCandidatePos(obj,detectParam, frames)
            %Method to perform localization on each plane for each frame
            %Check if some candidate exists already in the folder (previously saved)
            [run, candidate] = Core.MPLocMovie.existCandidate(obj.raw.movInfo.Path, '.mat');
            
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
            
        end
        
        function [candidate] = getCandidatePos(obj, frames)
            %Extract the position of the candidate of a given frame
            [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
            candidate = obj.candidatePos{idx};
            
            if isempty(candidate)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                
            end
        end
        
        function superResConsolidate(obj,roiSize,frames)
            %Consolidation refers to connect molecules that were localized
            %at similar position in different plane on a single frame.
            assert(~isempty(obj.calibrated),'Data should be calibrated to consolidate');
            assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before consolidation');
            
            %Check if some particles were saved already.
            [run, particle] = Core.MPLocMovie.existParticles(obj.raw.movInfo.Path, '.mat');
            
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
                    %#1 Extract Candidate Position for specific frame
                    [data] = obj.getFrame(idx);
                    [candidate] = obj.getCandidatePos(idx);
                    
                    if isempty(candidate)
                        
                        warning('Frame %d did not contain any candidate',idx);
                        particleList{idx} = [];
                        %particleList{idx}{1} = nan(5);
                        nParticles(idx) = 0;
                        
                    else
                        %#2 Consolidate the position of the given frame
                        %across plane
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
        
        function [traces,counter] = trackParticles(obj, trackParam)
            %track the particle in the Z direction (3rd dimension here)
            assert(~isempty(obj.calibrated),'Data should be calibrated to do ZzCalibrationration');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before zzCalibrationration');
            assert(~isempty(obj.particles), 'No particles found, please run superResConsolidate method before doing ZzCalibrationration');
            assert(isstruct(trackParam),'Tracking parameter is expected to be a struct with two field "euDistPx" and "ellip"')
            assert(and(isfield(trackParam,'euDistPx'),isfield(trackParam,'ellip')),...
                'Tracking parameter is expected to be a struct with two field "euDistPx" and "ellip"')
            %We copy the List as boolean to keep track of where there are
            %still particles left
            [listBool] = Core.MPLocMovie.copyList(obj.particles.List,1);
            %We copy as NaN for storage of the traces;
            [traces]   = Core.MPLocMovie.copyList(obj.particles.List,NaN);
            %We pick the first particle available
            [idx] = Core.MPLocMovie.pickParticle(listBool);
            counter = 1;
            errCount =1;
            while (idx)
                %loop until there is no particle (pickParticle return false)
                if errCount>1000
                    warning('While loop ran for unexpectedly longer time');
                    break;
                    
                end
                %Connect particles (cf consolidation but across frames
                [listIdx] = Core.MPLocMovie.connectParticles(obj.particles.List,listBool,idx, trackParam);
                %if the particle was connected in less than 5 frames we remove
                % its appearance from the list bool
                if length(listIdx) < 5
                    
                    [listBool] = Core.MPLocMovie.removeParticles(listBool,listIdx);
                    
                else
                    %Otherwise we store traces, increment counter and remove.
                    [traces]  = Core.MPLocMovie.storeTraces(traces,listIdx,counter);
                    counter = counter +1;
                    [listBool] = Core.MPLocMovie.removeParticles(listBool,listIdx);
                    
                end
                % We pick a new particle and start all over again
                [idx] = Core.MPLocMovie.pickParticle(listBool);
                errCount = errCount +1;
            end
            counter = counter -1;
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
            
            Thresh = 2*sqrt(2); %in Px As the superResCal is not
            %performed yet we allow 2px in both direction
            [checkRes1] = Core.MPLocMovie.checkEuDist(current(:,1:2),...
                next(:,1:2),Thresh);
            
            % Test ellipticity
            [checkRes2] = Core.MPLocMovie.checkEllipticity(current(:,3),...
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
        
        function [checkRes] = checkEllipticity(current, next, direction, thresh, eWeight)
            %Use to check if the Ellipticity make sense with what we
            %expect from the behavior of the PSFEngineering plate
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
        
        function [checkRes] = checkPlaneConfig(planeConfig,nPlanes)
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
        
        %method linked to consolidation along frames
        function listIdx    = connectParticles(List,listBool,idx, trackParam)
            %function connecting particles across frame by checking if two
            %particles from different frame are partners
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
                [checkRes] = Core.MPLocMovie.checkListBool(listBool,currentIdx(1)+1);
                nPartInFrame = length(checkRes);
                
                if ~all(checkRes==0)
                    %We use reshape to input the different particles of the next
                    %frame at once by storing them in the 3rd dimension
                    nextPart = List{currentIdx(1)+1};
                    nextPart = nextPart(logical(checkRes));
                    
                    shape = [size(nextPart{1},1) size(nextPart{1},2) nPartInFrame];%Unused particle are expected to be still 5x5
                    nextParts(:,:,1:nPartInFrame) = reshape(cell2mat(List{currentIdx(1)+1}'),shape);
                    nextParts = nextParts(:,:,logical(checkRes));
                    [isPart] = Core.MPLocMovie.isPartFrame(part2Track,nextParts,1,trackParam);
                    
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
                    
                    if counter == length(List)-1
                        
                        isPart = false;
                        
                    end
                    
                else
                    isPart = false;
                end
                listIdx(listIdx(:,1) == 0,:) = [];
                
            end
        end
        
        function [commonPlanes]  = findCommonPlanes(planeInCurrent,planeInNext)
            %Check how many common planes are between two candidates (one
            %of the criteria to determine if they can be partner)
            commonPlanes = zeros(size(planeInNext,1),2,size(planeInNext,2));
            
            for i = 1 : size(planeInNext,2)
                
                commonPlanes(:,1,i) = ismember(planeInCurrent,planeInNext(:,i));
                commonPlanes(:,2,i) = ismember(planeInNext(:,i),planeInCurrent);
                
            end
            commonPlanes = logical(squeeze(commonPlanes));
        end
        
        function [isPart]   = isPartFrame(current, next, direction, trackParam)
            %This function is designed to have PSFE plate ON
            assert(abs(direction) == 1, 'direction is supposed to be either 1 (up) or -1 (down)');
            assert(size(current,2) == size(next,2), 'Dimension mismatch between the tail and partner to track');
            
            if ~all(isnan(next))
                %ZStack, consolidation between frame
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
                    commonPlanes = Core.MPLocMovie.findCommonPlanes(current(:,end),squeeze(next(:,end,:)));
                    roughcheck2  = sum(commonPlanes,1)>=nConsistentPlanes;
                    if all(roughcheck2 ==0)
                        disp('Less than 2 planes in common, breaking out');
                    else
                        
                        for i = 1 : size(next,3)
                            % Test Euclidian distance
                            Thresh = trackParam.euDistPx; %in px
                            [checkRes1] = Core.MPLocMovie.checkEuDist(current(commonPlanes(:,1,i),1:2),...
                                squeeze(next(commonPlanes(:,2,i),1:2,i)),Thresh);
                            
                            % Test ellipticity
                            eWeight = [1 2 3 2 1];
                            thresh = trackParam.ellip;
                            [checkRes2] = Core.MPLocMovie.checkEllipticity(current(commonPlanes(:,1,i),3),...
                                squeeze(next(commonPlanes(:,2,i),3,i)),direction,thresh,eWeight(commonPlanes(:,1,i)));
                            %To be a particle, we want the position to be
                            %consistent in at least 2 planes and
                            %ellipticity to pass the test.
                            isPart(i) = and(length(find(checkRes1))>=nConsistentPlanes, checkRes2);
                            
                        end
                        
                    end
                end
                
                isPart = logical(isPart);
                
            else
                isPart = false(size(next,1),1);
            end
        end
        
        function [listCopy] = copyList(List, filling)
            %Small function that copy a given list (cell array) filling it
            %with the input filling (e.g. NaN, 0, 1,...)
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
        
        function [checkRes] = checkListBool(listBool, idx)
            %Check if there are still candidate in the list that were not
            %yet connected or removed.
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
        
        function [idx]      = pickParticle(listBool)
            %This function will pick a particle and return false if there
            %is none. It goes through the given list and search for the
            %first non-false idx and returns it.
            
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
        
        function [traces]   = storeTraces(traces,listIdx,numPart)
            %Store the list of idx as a traces (which has the same format
            %as the particle list
            for i = 1 : size(listIdx,1)
                
                traces{listIdx(i,1)}{listIdx(i,2)} = numPart;
                
            end
            
        end
        
        function [listBool] = removeParticles(listBool,listIdx)
            %Remove the particle by putting false in the boolean copy in
            %the given indices of the indices list (listIdx)
            for i = 1 : size(listIdx,1)
                
                listBool{listIdx(i,1)}{listIdx(i,2)} = false;
                
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
                
                position = zeros(200,3);
                [volIm] = obj.getFrame(frames(i));
                nameFields = fieldnames(volIm);
                
                for j = 1:length(nameFields)
                    %localization occurs here
                    [ pos, ~, ~ ] = Localization.smDetection( double(volIm.(nameFields{j})),...
                        delta, FWHM_pix, chi2 );
                    startIdx = find(position==0,1,'First');
                    pos(:,3) = j;
                    position(startIdx:startIdx+size(pos,1)-1,:) = pos;
                    
                end
                
                idx = find(position==0,1,'First');
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
        
        %Methods linked to consolidation
        function [finalCandidate] = consolidatePos(obj, data, frameCandidate, roiSize)
            %SuperRes fitting and GLRT focus metric are determined here and
            %fused in candMetric
            [candMetric] = obj.superResLocFit(data,frameCandidate,roiSize);
            %Calculate a focus metric (FM) combining ellipticity and GLRT FM.
            [corrEllip, focusMetric] = Localization.calcFocusMetric(candMetric(:,3),candMetric(:,4));
            
            %reformating to keep the same format as how the data is saved
            %later
            candMetric = [candMetric candMetric(:,end)];
            candMetric(:,6) = candMetric(:,5);
            candMetric(:,5) = focusMetric;
            candMetric(:,4) = [];
            focusMetric((1-corrEllip)>0.3) = NaN;
            %Plane Consolidation occur here
            [finalCandidate] = obj.planeConsolidation(candMetric,focusMetric);
            
            %we delete empty cells from the array
            idx2Empty = cellfun(@isempty,finalCandidate);
            finalCandidate(idx2Empty(:,1),:) = [];
            
        end
        
        function [candMet] = superResLocFit(obj,data,frameCandidate,roiSize)
            %Candidate metric are determined here (x,y,e,+focusmetric)
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
                        %LRT focus metric
                        [LRT,~] = Localization.likelihoodRatioTest(ROI,sig,[row col]);
                        rowPos = planeCandidate(k,1) + row;
                        colPos = planeCandidate(k,2) + col;
                        
                        [grad,~] = imgradient(ROI);
                        Grad = max(max(grad,[],2),[],1);
                        cM(k,:) = [rowPos colPos e LRT Grad];
                        
                    end
                    %candidate Metric: phasor localization Res, Likelihood Res,
                    %indexCandidate in plane, Plane.
                    candMet(currentk:currentk+k-1,:) = [cM j*ones(k,1)];
                    %candidateMetric: [row col e LRT Grad J]
                    currentk = currentk+k;
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
                currentPlane = candMet(idx,end);
                
                %Check which planes are to be checked (currently 2 planes
                %above and 2 planes below the given plane
                planes2Check = currentPlane-2:currentPlane-1;
                planes2Check = planes2Check(planes2Check>0);
                planes2Check = [planes2Check currentPlane+1:currentPlane+2];
                planes2Check = planes2Check(planes2Check<nPlanes+1);
                currentCand = candMet(idx,:);
                direction = 1;%Start by checking above
                
                particle = nan(5,6);
                particle(3,:) = currentCand;
                nCheck = length(planes2Check);
                for i = 1:nCheck
                    
                    cand = candMet(candMet(:,end) == planes2Check(i),:);
                    if(planes2Check(i) > currentPlane)
                        direction = -1;%check below
                    end
                    
                    [isPart] = Core.MPLocMovie.isPartPlane(currentCand,cand,direction);
                    if ~all(isPart ==0)
                        id = cand(isPart,end)-currentCand(end);
                        particle(round(length(currentCand)/2)+id,:) = cand(isPart,:);
                    end
                    
                end
                
                %Check if the resulting configuration of the plane make
                %sense e.g. no hole in the configuration
                planeConfig = particle(:,end);
                [checkRes] = Core.MPLocMovie.checkPlaneConfig(planeConfig,nPlanes);
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
end