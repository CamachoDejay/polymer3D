classdef MPParticleMovie < Core.MPMovie
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'protected')
        
        candidatePos
        unCorrLocPos
        corrLocPos
        particles
        intData
        
    end
    
    methods
        function obj = MPParticleMovie(raw,cal,info)
            
            obj  = obj@Core.MPMovie(raw,cal,info);
            
        end
        
        function set.candidatePos(obj,candidatePos)
            
            assert(iscell(candidatePos), 'Expecting a cell containing candidate positions');
            obj.candidatePos = candidatePos;
            
        end
        
        function findCandidatePos(obj,detectParam, frames)
            %Method to perform localization on each plane for each frame
            %Check if some candidate exists already in the folder (previously saved)
            [run, candidate] = obj.existCandidate(obj.raw.movInfo.Path, '.mat');
            
            if run
                switch nargin
                    case 2
                        
                        frames = 1: obj.calibrated.nFrames;
                        disp('Running detection on every frame');
                        
                    case 3
                        
                        [frames] = obj.checkFrame(frames,obj.calibrated.nFrames);
                        
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
        
        function consolidatePlanes(obj,roiSize,frames,consThresh)
            %Consolidation refers to connect molecules that were localized
            %at similar position in different plane on a single frame.
            assert(~isempty(obj.calibrated),'Data should be calibrated to consolidate');
            assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before consolidation');
            assert(~isempty(obj.unCorrLocPos),'Localization needs to be performed before consolidation');
           
            %Check if some particles were saved already.
            [run, particle] = obj.existParticles(obj.raw.movInfo.Path, '.mat');
            
            if run
                %Check the number of function input
                switch nargin
                    case 1
                        
                        frames = 1: obj.calibrated.nFrames;
                        roiSize = 6;
                        disp('Running consolidation on every frame with roi of 6 pixel');
                        consThresh = 3;
                    case 2
                        
                        frames = 1: obj.calibrated.nFrames;
                        disp('Running consolidation on every frame')
                        consThresh = 3;
                    case 3
                        
                        [frames] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                        assert(min(size(roiSize))==1,'RoiSize is expected to be a single number')
                        assert (isnumeric(roiSize),'RoiSize is expected to be a single number');
                        consThresh = 3;
                    case 4
                        [frames] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
                        assert(min(size(roiSize))==1,'RoiSize is expected to be a single number')
                        assert (isnumeric(roiSize),'RoiSize is expected to be a single number');
                        
                    otherwise
                        
                        error('Something wrong with number of input');
                        
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
                          switch obj.info.zMethod
                              case 'Intensity'
                                 %we do not do anythin at the moment.
                                 focusMetric = abs(fCandMet.fMetric);
                              case 'PSFE' 
                                [corrEllip, focusMetric] = Localization.calcFocusMetric(fCandMet.ellip,fCandMet.fMetric);
                          end
                            %reformating to keep the same format as how the data is saved
                            %later
                            fCandMet.fMetric = focusMetric;
                            
                            %focusMetric((1-corrEllip)>0.3) = NaN;
                            
                            %Plane Consolidation occur here
                            [part] = obj.planeConsolidation(fCandMet,focusMetric,consThresh);

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
            [idx] = obj.checkFrame(frames,obj.calibrated.nFrames(1));
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
            
            nImages = size(frame,3);
           
            nsFig = 2;
            
            candidate = obj.getCandidatePos(idx);
            rowPos    = candidate.row;
            colPos    = candidate.col;
            planeIdx  = candidate.plane;
            
            h = figure(2);
            h.Name = sprintf('Frame %d',idx);
            for i = 1:nImages
                
                subplot(2,nImages/nsFig,i)
                hold on
                imagesc(frame(:,:,i))
                plot(colPos(planeIdx==i),rowPos(planeIdx==i),'g+')
                axis image;
                grid on;
                a = gca;
                a.XTickLabel = [];
                a.YTickLabel = [];
                a.GridColor = [1 1 1];
                title({['Plane ' num2str(i)],sprintf(' Zpos = %0.3f',obj.calibrated.oRelZPos(i))});
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
                            if(~isempty(currPart(currPart.plane == i,:)))
                                part2Plot = currPart(currPart.plane == i,:);
                                plot(part2Plot.col,part2Plot.row,'o',...
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
                        idx2NaN = isnan(currPart.row);
                        currPart(idx2NaN,:) = [];
                        planes = currPart.plane;
                        figure(20+i)
                        hold on
                        for j = 1 : length(planes)
                            jdx = planes(j);
                            currFrame = frame.(sprintf('plane%d',jdx));
                            ROI = EmitterSim.getROI(currPart.col(j), currPart.row(j),...
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
        
        function [candidateList] = planeConsolidation(obj,candMet,focusMetric,consThresh)
            %Loop through all candidate of a given frame and match them
            %between frame until none can be match or all are matched.
            nPlanes = obj.calibrated.nPlanes;
            counter = 1;
            nPart = 0;
            maxIt = size(candMet,1);
            zMethod = obj.info.zMethod;
            candidateList = cell(max(size(find(~isnan(focusMetric)))),1);
            %continue until the list is not empty
            while and(~isempty(focusMetric), ~isnan(nanmax(focusMetric)))
                
                if counter> maxIt
                    
                    error('While loop ran for an unexpectedly long time, something might be wrong');
                    
                end
                
                %Find candidate in best focus
                [~,idx] = max(focusMetric);
                currentPlane = candMet.plane(idx);
                
                switch nPlanes
                    case 1
                        planes2Check = [];
                        
                    otherwise
                        
                        %Check which planes are to be checked (currently 2 planes
                        %above and 2 planes below the given plane
                        planes2Check = currentPlane-nPlanes:currentPlane-1;
                        planes2Check = planes2Check(planes2Check>0);
                        planes2Check = [planes2Check currentPlane+1:currentPlane+nPlanes];
                        planes2Check = planes2Check(planes2Check<nPlanes+1);
                        
                end
                currentCand = candMet(idx,:);
                direction = -1;%Start by checking above
                
                particle = array2table(nan(nPlanes,size(currentCand,2)));
                particle.Properties.VariableNames = currentCand.Properties.VariableNames;
                
                particle(currentCand.plane,:) = currentCand;
                nCheck = length(planes2Check);
                camConfig = obj.calibrated.camConfig;
                for i = 1:nCheck
                    
                    cand = candMet(candMet.plane == planes2Check(i),:);
                    if(planes2Check(i) > currentPlane)
                        direction = +1;%check below (Plane 1 is the uppest plane 8 is lowest)
                    end
                    
                    [isPart] = Core.MPParticleMovie.isPartPlane(currentCand,cand,direction,consThresh,zMethod);
                    if ~all(isPart ==0)
                        id = cand.plane(isPart);
                        particle(id,:) = cand(isPart,:);
                    end
                    
                end
               
                %We remove the particle(s) from the list
                focusMetric(ismember(candMet(:,1), particle(:,1))) = [];
                candMet(ismember(candMet(:,1), particle(:,1)),:) = [];
                %format particle to be the same as before:
                [particle] = obj.makeParticle(particle);
                
                 %Check if the resulting configuration of the plane make
                %sense e.g. no hole in the configuration
                 planeConfig = particle.plane;
                [checkRes] = Core.MPParticleMovie.checkPlaneConfig(planeConfig,nPlanes,camConfig,zMethod);
                
                %Store
                if checkRes
                    
                    nPart = nPart +1;
                    %store particle in a new list
                    candidateList{nPart,1} = particle;
                    
                else
                    %Otherwise we remove it from the best focus search list
                    %by putting focus metric to NaN
                    %focusMetric(idx) = NaN;
                    
                end
                
                counter = counter+1;
                
            end
        end
                
    end
     methods (Static)
       
        function [isPart]   = isPartPlane(current, next, direction,consThresh,zMethod)
            %This function aim at determining whether a candidate from one
            %plane and the another are actually the same candidate on
            %different plane or different candidate. The decision is based
            %on threshold on localization distance, ellipticity and focus
            %metric.
            
            %This function is designed to have PSFE plate ON
            assert(abs(direction) == 1, 'direction is supposed to be either 1 (up) or -1 (down)');
            assert(size(current,2) == size(next,2), 'Dimension mismatch between the tail and partner to track');
            
            thresh = consThresh;
            [checkRes1] = Core.MPParticleMovie.checkEuDist([current.row, current.col],...
                [next.row, next.col],thresh);
            
            if strcmp(zMethod,'PSFE')
             % Test ellipticity
                [checkRes2] = Core.MPParticleMovie.checkEllipticity(current.ellip,...
                next.ellip,direction);
            
            elseif strcmp(zMethod,'Intensity')
                
                checkRes2 = checkRes1;
                
            else
                
                error('Unknown Z method for consolidation');
                
            end
            
            % Test focus Metric
            maxExpFM = current.fMetric+0.1*current.fMetric;
            checkRes3 = next.fMetric < maxExpFM;
            
            %isPart will only be true for particle that passes the 3 tests
            isPart = checkRes1.*checkRes2.*checkRes3;
            
            if(length(find(isPart))>1)
                
                warning('Could not choose which particle was the partner of the requested particle, killed them both');
                isPart(isPart==1) = 0;
            end
            
            if isempty(isPart)
                isPart = false;
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
        
        function [checkRes] = checkPlaneConfig(planeConfig,nPlanes,camConfig,zMethod)
            %Here we will check that the consolidation found based on the
            %best focused particle make sense with what we would expect and
            %also that we have enough planes.
            if strcmp(zMethod,'Intensity')
                testPlanes = sum(~isnan(planeConfig))>=2;
            elseif strcmp(zMethod,'PSFE')
                switch nPlanes
                    case 1
                        nPlanesEdgeFrange = 1;
                        nPlanesFullRange = 1;
                        nPlanesInterleaved = 1;
                    case 2

                        nPlanesEdgeFrange = 1;
                        nPlanesFullRange = 1;
                        nPlanesInterleaved = 1;

                    case 4

                        nPlanesEdgeFrange = 1;
                        nPlanesFullRange = 2;
                        nPlanesInterleaved = 2;

                    case 8

                        nPlanesEdgeFrange = 1;
                        nPlanesEdgeInterleaved = 2;
                        nPlanesFullRange = 2;
                        nPlanesInterleaved = 3;
                    otherwise
                        error('Unknown number of planes, only expect 1,2,4,8')
                end
                %Let us test that we have consolidate the particle in at least
                %3 Planes
                isEdgePlane = or(~isempty(find(planeConfig==1,1)),~isempty(find(planeConfig==8,1)));


                switch camConfig
                    case 'fullRange'
                        if isEdgePlane

                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesEdgeFrange;

                        else
                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesFullRange;
                        end

                    case 'interleaved'
                        if isEdgePlane
                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesEdgeInterleaved;
                        else

                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesInterleaved;
                        end
                    case 'equal'
                        if isEdgePlane

                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesEdgeFrange;

                        else
                            testPlanes = length(find(~isnan(planeConfig)==true)) >= nPlanesFullRange;
                        end

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
        
        function [int,SNR]  = getIntensity(ROI,sig)
            %extract central position
            center = [ceil(size(ROI,1)/2),ceil(size(ROI,2)/2)];
            rowPos = center(1);
            colPos = center(2);
            %we integrate at 3*sigma (take ROI
            roiSignal = ceil(sig);
            
            %get the idx for the ROI to integrate
            rowIdx = rowPos-roiSignal(1):rowPos+roiSignal(1);
            colIdx = colPos-roiSignal(2):colPos+roiSignal(2);
            
            %Pixel to integrate for signal
            px2SumInt = ROI(rowIdx,colIdx);
            
             %get background pixels
            bkg = ROI;
            bkg(rowIdx,colIdx) = 0;
            px2SumBkg = bkg(bkg~=0);
            bkg = mean(px2SumBkg);
            bkgVar = std(px2SumBkg);
            %calculate signal
            int = px2SumInt - bkg;
            int = sum(sum(int));
            %SNR = max(max(px2SumInt))/bkgVar;
            SNR = sqrt(int);
            if or(SNR<0,int<0)
                int = sum(sum(px2SumInt));
                SNR = sqrt(int);
            end
           
        end
       
     end
     
     methods (Access = protected)
        %method linked to candidate
        function [run,candidate] = existCandidate(obj,Path,ext)
            
            runMethod = obj.info.runMethod;
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'candidatePos')==true)
                        candidate = load([file2Analyze(1).folder filesep 'candidatePos.mat']);
                        candidate = candidate.candidate;
                        run = false;
                    else
                
                        run = true;
                        candidate =[];
                
                    end
                case 'run'
                    run = true;
                    candidate =[];
            end
        end
        
        %method linked to fitting
        function [run,SRLocPos] = existLocPos(obj,Path,ext) 
            runMethod = obj.info.runMethod;
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'SRLocPos')==true)
                        SRLocPos = load([file2Analyze(1).folder filesep 'SRLocPos.mat']);
                        name = fieldnames(SRLocPos);
                        SRLocPos = SRLocPos.(name{1});
                        run = false;
                    else
                
                         run = true;
                        SRLocPos =[];
                
                    end
                case 'run'
                     run = true;
                     SRLocPos =[];
            end    
        end
        
        %method Linked to particles/planeConsolidation
        function [run, particle] = existParticles(obj,Path, ext)
            
            runMethod = obj.info.runMethod;
            switch runMethod
                case 'load'
                    [file2Analyze] = Core.Movie.getFileInPath(Path, ext);
                    %Check if some candidate were already stored
                    if any(contains({file2Analyze.name},'particle')==true)
                        particle = load([file2Analyze(1).folder filesep 'particle.mat']);
                        particle = particle.particle;
                        run = false;
                    else
                
                        run = true;
                        particle = [];
                
                    end
                    
                case 'run'
                    
                     run = true;
                     particle = [];
                     
            end
        end
        
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
                nPlanes = size(volIm,3);
                
                for j = 1:nPlanes
                    currentIM = volIm(:,:,j);
                    %localization occurs here
                    [ pos, meanFAR, ~ ] = Localization.smDetection(currentIM,...
                        delta, FWHM_pix, chi2 );
                    if ~isempty(pos)
                        startIdx = find(position.row==0,1,'First');
                        if isempty(startIdx)
                            startIdx = length(position.row)+1;
                        end
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
                sigSetup = [obj.info.sigma_px obj.info.sigma_px];
                
            for i = 1:size(frameCandidate,1)
                
                plane = frameCandidate.plane(i);
                planeData = double(data(:,:,plane));
                %Get the ROI
                [roi_lims] = EmitterSim.getROI(frameCandidate.col(i), frameCandidate.row(i),...
                    delta, size(planeData,2), size(planeData,1));
                ROI = planeData(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2));
                roiSize = size(ROI);
                if roiSize(1)==roiSize(2)
                    if strcmpi(obj.info.fitMethod,'phasor')
                        %Phasor fitting to get x,y,e
                        [row,col,e,magX,magY] = Localization.phasor(ROI);
                        rowPos = round(frameCandidate.row(i)) + row;
                        colPos = round(frameCandidate.col(i)) + col;

                    elseif strcmpi(obj.info.fitMethod,'Gauss')
                        [X,Y] = meshgrid(frameCandidate.col(i)-delta:frameCandidate.col(i)+...
                            delta,frameCandidate.row(i)-delta:frameCandidate.row(i)+delta);
                        domain(:,:,1) = X;
                        domain(:,:,2) = Y;

                        %Gauss (slower)
                        [gPar] = Localization.Gauss.MultipleFitting(ROI,frameCandidate.col(i),...
                            frameCandidate.row(i),domain,1);%data,x0,y0,domain,nbOfFit
                        colPos = gPar(5); %Should be directly the position of the particle as we
                            %gave above the domain of the ROI in the space of the image
                        rowPos = gPar(6);

                        row = rowPos - round(frameCandidate.row(i));
                        col = colPos - round(frameCandidate.col(i));

                        e = gPar(3)/gPar(2);

                        magX = 0;
                        magY = 0;
                    else
                        %Phasor fitting to get x,y,e
                        [row,col,e,magX,magY] = Localization.phasor(ROI);
                        rowPos = round(frameCandidate.row(i)) + row;
                        colPos = round(frameCandidate.col(i)) + col;
                    end

                     %LRT focus metric
                    [fMetric,~] = Localization.likelihoodRatioTest(ROI,sigSetup,[row col]);

                    if magX>=magY
                        sig(1) = sigSetup(1) * magX/magY;
                        sig(2) = sigSetup(2);
                    else
                        sig(1) = sigSetup(1);
                        sig(2) = sigSetup(2) * magY/magX;
                    end

                    [int,SNR] = obj.getIntensity(ROI,sigSetup);
                    %LRT focus metric
                    [gFitMet,~] = Localization.likelihoodRatioTest(ROI,sig,[row col]);

                    %storing info
                    candMet.row(i) = rowPos;
                    candMet.col(i) = colPos;
                    candMet.z(i) = 0;
                    candMet.ellip(i) = e;
                    candMet.magX(i) = magX;
                    candMet.magY(i) = magY;
                    candMet.intensity(i) = int;
                    candMet.SNR(i) = SNR;
                    candMet.meanFAR(i) = frameCandidate.meanFAR(i);
                    candMet.fMetric(i) = fMetric;
                    candMet.gFitMet(i) = gFitMet;
                    candMet.plane(i) = plane;
                else
                    
                    %storing info
                    candMet.row(i)          = NaN;
                    candMet.col(i)          = NaN;
                    candMet.z(i)            = NaN;
                    candMet.ellip(i)        = NaN;
                    candMet.magX(i)         = NaN;
                    candMet.magY(i)         = NaN;
                    candMet.intensity(i)    = NaN;
                    candMet.SNR(i)          = NaN;
                    candMet.meanFAR(i)      = NaN;
                    candMet.fMetric(i)      = NaN;
                    candMet.gFitMet(i)      = NaN;
                    candMet.plane(i)        = NaN;
                end
            end
            %remove NaN's  
            idx = isnan(candMet.row);
            candMet(idx,:) = [];

        end
        
        function [doAvg]  = checkDoAverage(obj,ellip)
            camConfig = obj.calibrated.camConfig;
            switch camConfig
                case 'fullRange'
                    if and(ellip>0.9,ellip<1.1)
                        doAvg = false;
                    else
                        doAvg = true;
                    end

                case 'interleaved'

                        doAvg = true;

                case 'equal'
                    if and(ellip > 0.9, ellip < 1.1)
                        
                        doAvg = false;
                    
                    else
                        
                        doAvg = true;
                    
                    end

                otherwise
                    error('unknown camera config');
            end
        end
        
        function [newPart] = makeParticle(~,particleData)
            newPart = array2table(nan(5,size(particleData,2)));
            newPart.Properties.VariableNames = particleData.Properties.VariableNames;
            %store best focus in center
            [~,idx] = nanmax(particleData.fMetric);
            if idx-2 > 0
                newPart(1,:) = particleData(idx-2,:);
            end
            
            if idx-1 > 0
                newPart(2,:) = particleData(idx-1,:);
            end
            
            newPart(3,:) = particleData(idx,:);
            
            if idx+1 <= 8
                newPart(4,:) = particleData(idx+1,:);
            end
            
            if idx+2 <= 8
                newPart(5,:) = particleData(idx+2,:);
            end
        
        end
       
     end
end

