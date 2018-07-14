classdef mpLocMovie < Core.mpMovie
    %mpLocMovie hold all the method and info necessary for localization,
    %expension of this object will be first centered around Tracking but
    %it will also be easy to extend to STORM/PALM.
    
    properties
        
        particles
        candidatePos
        
    end
    
    methods
        
        function obj = mpLocMovie(raw, cal)
            
            obj  = obj@Core.mpMovie(raw,cal);
         
        end
        
        function set.candidatePos(obj,candidatePos)
            
            assert(iscell(candidatePos), 'Expecting a cell containing candidate positions');
            obj.candidatePos = candidatePos;
            
        end
            
        function [particle] = getParticles(obj,frames)
          %GetParticles
            [idx] = obj.checkFrame(frames);
            particle = obj.particles.List{idx};
            
            if isempty(particle)
                
                warning('There was no particle found in this frames, check than you ran superResConsolidate method on this frame beforehand');
            
            end
        end
        
        function findCandidatePos(obj,detectParam, frames)
            %Method to perform localization on each plane for each frame
            assert(~isempty(obj.info), 'Missing information about setup to be able to find candidates, please use giveInfo method first');
            assert(nargin>1,'not enough input argument');
            [file2Analyze] = getFileInPath(obj, obj.raw.movInfo.Path, '.mat');
            %Check if some candidate exists already (previously saved)
            [run, candidate] = obj.existCandidate(file2Analyze);
            
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
            [idx] = obj.checkFrame(frames);
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
            
            [file2Analyze] = getFileInPath(obj, obj.raw.movInfo.Path, '.mat');            
            [run, particle] = obj.existParticles(file2Analyze);
            
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
                %Consolidation occurs here
                for i = 1 : 1:nFrames
                    disp(['Consolidating frame ' num2str(i) ' / ' num2str(nFrames)]);
                    idx = frames(i);
                    %#1 Extract Candidate Position for specific frame
                    [data] = obj.getFrame(idx);
                    [candidate] = obj.getCandidatePos(idx);

                    if isempty(candidate)

                        warning('Frame %d did not contain any candidate',idx);
                        particleList{idx} = nan(5);
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
        
        function showCandidate(obj,idx)
            %Display Candidate
            assert(length(idx)==1, 'Only one frame can be displayed at once');
            [idx] = obj.checkFrame(idx);
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
        
        function [checkRes] = checkEuDist(obj,current,next,Thresh)
            %Use to check if the Euclidian distance is within reasonable
            %range
            EuDist = sqrt((current(:,1) - next(:,1)).^2 +...
                            (current(:,2) - next(:,2)).^2);
            checkRes = EuDist < Thresh;
            
            if isempty(checkRes)
                checkRes = false;
            end
            
        end
        
        function [checkRes] = checkEllipticity(obj, current, next, direction, thresh, eWeight)
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
        
        function [run, particle] = existParticles(obj,file2Analyze)
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
                end
            else
                run = true;
                particle = [];
            end    
            
        end
        
        function [commonPlanes] = findCommonPlanes(~,planeInCurrent,planeInNext)
            %Check how many common planes are between two candidates (one
            %of the criteria to determine if they can be partner)
            commonPlanes = zeros(size(planeInNext,1),2,size(planeInNext,2));
            
            for i = 1 : size(planeInNext,2)
                
                commonPlanes(:,1,i) = ismember(planeInCurrent,planeInNext(:,i));
                commonPlanes(:,2,i) = ismember(planeInNext(:,i),planeInCurrent);
                
            end
            commonPlanes = logical(squeeze(commonPlanes));
        end
        
    end
    
    methods (Access = private)
        
        %Methods linked to Candidate
        function [run,candidate] = existCandidate(obj,file2Analyze)
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
                
            end
            
            else
                
                run = true;
                candidate =[];
            end    
        end
        
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
                candidate{frames(i)} = position(1:idx-1,:);
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
        
        function [isPart]   = isPartPlane(obj, current, next, direction)
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

            isPart = logical(isPart);         
                
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
                    
                    [isPart] = obj.isPartPlane(currentCand,cand,direction);
                    if ~all(isPart ==0)
                        id = cand(isPart,end)-currentCand(end);
                        particle(round(length(currentCand)/2)+id,:) = cand(isPart,:);
                    end
                    
                end

                %Check if the resulting configuration of the plane make
                %sense e.g. no hole in the configuration
                planeConfig = particle(:,end);
                [checkRes] = obj.checkPlaneConfig(planeConfig);
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