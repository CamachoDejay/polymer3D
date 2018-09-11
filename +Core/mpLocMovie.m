classdef MPLocMovie < Core.MPParticleMovie
    %MPLocMovie hold all the method and info necessary for localization,
    %expension of this object will be first centered around Tracking but
    %it will also be easy to extend to STORM/PALM.
    
    properties (SetAccess = 'protected')
        
        SRCal
        ZCal
        corrLocPos
    
    end
    
    methods
        
        function obj = MPLocMovie(raw, MPCal, SRCal, zCal)
            
            obj  = obj@Core.MPParticleMovie(raw,MPCal);
            
            switch nargin
                
                case 1
                    error('MPCal is required to create MPLocMovie')
                case 2
                    obj.SRCal = [];
                    obj.ZCal  = [];
                case 3
                    obj.SRCal = SRCal;
                    obj.ZCal = [];
                case 4
                    obj.SRCal = SRCal;
                    obj.ZCal = zCal;
                otherwise
                    error('Too many input arguments');
            end
        end
        
        function set.SRCal(obj,SRCal)
            assert(isfolder(SRCal), 'The given path is not a folder');
            
            %Check Given path
            [file2Analyze] = Core.Movie.getFileInPath(SRCal,'SRCalibration.mat');
            
            if isempty(file2Analyze)
                error('No SR calibration file found in the given folder');
            else
                fileName = [file2Analyze.folder filesep file2Analyze.name];
                cal = load(fileName);
                field = fieldnames(cal);
                cal = cal.(field{1});
                assert(and(isfield(cal,'translation'),isfield(cal,'rotation')),...
                    'SR calibration is supposed to have a tranlsation and rotation field');
                
                obj.SRCal = cal; 
            end
            
        end
        
        function set.ZCal(obj,zCal)
            
            assert(isfolder(zCal), 'The given path is not a folder');
            
            %Check Given path
            [file2Analyze] = Core.Movie.getFileInPath(zCal,'zCalibration.mat');
            
            if isempty(file2Analyze)
                error('No z calibration file found in the given folder');
            else
                fileName = [file2Analyze.folder filesep file2Analyze.name];
                cal = load(fileName);
                field = fieldnames(cal);
                cal = cal.(field{1});
                assert(isstruct(cal),'zCalibration is supposed to be in cells format');
                assert(and(isfield(cal,'fitZParam'),isfield(cal,'calib')),...
                    'Something is wrong in the fields of your Z calibration');
                
                obj.ZCal = cal; 
            end
        end
        
        function applyCorr(obj,rot,refPlane)
            
            %apply SRCal
            obj.applySRCal(rot,refPlane);
            
            %transform ellipticity into Z
            obj.applyZCal;
            
        end
        
        function applySRCal(obj, rot, refPlane)
            assert(~isempty(obj.localizedPos),'You need to find candidate and SR Localized them before applying corrections');
            assert(~isempty(obj.SRCal),'SR Calibration needed to correct the data');
            
            if nargin <2
                refPlane = 5;
            end
            if isempty(obj.corrLocPos)
                obj.corrLocPos = obj.localizedPos;
            end
            
            data = obj.localizedPos;
            rotMat = obj.SRCal.rotation;
            transMat = obj.SRCal.translation;
           % correctedData = cell(size(data));
            disp(['Applying SR calibration...']);
            for i = 1 : length(data)
                currData = data{i};
                currPlanes = unique(currData(:,end));
                for j = 1 : length(currPlanes)
                    currentPlane = currPlanes(j);
                    data2Corr = currData(currData(:,end)==currentPlane,1:2);
                    
                    [corrData] = obj.applyTrans(data2Corr,transMat,refPlane,currentPlane);
                    
                    if rot 
                        [corrData] = obj.applyRot(corrData, rotMat,refPlane,currentPlane);
                    end
                    %we store the corrected data
                    obj.corrLocPos{i}(currData(:,end)==currentPlane,1:2) = corrData;
                    
                end
                %correctedData{i}(:,3) = data{i}(:,3);
            end
            disp('========> DONE ! <=========');
           % obj.corrLocPos = correctedData;
            
        end
        
        function applyZCal(obj)
            
            assert(~isempty(obj.localizedPos),'You need to find candidate and SR Localized them before applying corrections');
            assert(~isempty(obj.ZCal),'Z Calibration needed to correct the data');
            
            if isempty(obj.corrLocPos)
                obj.corrLocPos = obj.localizedPos;
                warning('Z calibration is currently being applied on non-SRCorrected (X-Y) data');
            end
            
            data = obj.localizedPos; 
            zCal = obj.ZCal;
            
            %Here we translate ellipticity into z position based on
            %calibration
            disp('Applying Z Calibration');
            for i = 1 : length(data)
                currData = data{i};
                nPos = size(currData,1);
                
                for j = 1 : nPos
                    
                    currentEllip = currData(j,3);
                    currentPlane = currData(j,end);
                    [zPos] = getZPosition(obj,currentEllip,zCal,currentPlane);
                    obj.corrLocPos{i}(j,3) = zPos;
                    
                end
               
            end
            
            %Here we translate the ellipticity range into zRange for each
            %plane
            
            ellipRange = zCal.fitZParam.ellipRange;
            nPlanes = obj.calibrated.nPlanes;
            zRange = cell(nPlanes,1);
            for i = 1 : nPlanes
                zRange{i} = obj.getZRange(ellipRange,zCal,i);
            end
            disp('=======> DONE ! <========');
        end
        
        function [locPos] = getCorrLocPos(obj,frames)
             %Extract the position of the candidate of a given frame
            [idx] = Core.Movie.checkFrame(frames,obj.raw.maxFrame(1));
            locPos = obj.corrLocPos{idx};
            
            if isempty(locPos)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
                
            end
        end
           
%         function consolidatePlanes(obj,frames)
%             
%             %Consolidation refers to connect molecules that were localized
%             %at similar position in different plane on a single frame.
%             assert(~isempty(obj.calibrated),'Data should be calibrated to consolidate');
%             assert(~isempty(obj.info),'Information about the setup are missing to consolidate, please fill them in using giveInfo method');
%             assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before consolidation');
%             assert(~isempty(obj.localizedPos),'Localization needs to be performed before consolidation');
%            
%             %Check if some particles were saved already.
%             [run, particle] = Core.MPParticleMovie.existParticles(obj.raw.movInfo.Path, '.mat');
%             
%             if run
%                 %Check the number of function input
%                 switch nargin
%                     case 1
%                         
%                         frames = 1: obj.calibrated.nFrames;
%                         
%                         disp('Running consolidation on every frame...');
%                         
%                     case 2
%                         
%                        [frames] = Movie.checkFrame(frames,obj.raw.maxFrame(1));
%                        
%                     otherwise
%                         
%                         error('Something wrong with input');
%                         
%                 end
%                 
%                 nFrames = length(frames);
%                 %allocate for storage
%                 particleList = cell(1,obj.raw.maxFrame(1));
%                 nParticles = zeros(1,obj.raw.maxFrame(1));
%                 idx2TP = zeros(1,obj.raw.maxFrame(1));
%                 h = waitbar(0,'Consolidating candidate ...');
%                 
%                 %Consolidation occurs here
%                 for i = 1 : 1:nFrames
%                     disp(['Consolidating frame ' num2str(i) ' / ' num2str(nFrames)]);
%                     idx = frames(i);
%                     %#1 Extract localized Position for specific frame
%                     [fCandMet] = obj.getLocPos(idx);
%                     [candMet]  = obj.getCorrLocPos(idx);
%                     if isempty(fCandMet)
%                         
%                         warning('Frame %d did not contain any localized positions',idx);
%                         particleList{idx} = [];
%                         %particleList{idx}{1} = nan(5);
%                         nParticles(idx) = 0;
%                         
%                     else
%                         %#2 Consolidate the position of the given frame
%                         %across plane
%                           %Calculate a focus metric (FM) combining ellipticity and GLRT FM.
%                             [corrEllip, focusMetric] = Localization.calcFocusMetric(fCandMet(:,3),fCandMet(:,4));
% 
%                             %reformating to keep the same format as how the data is saved
%                             %later
%                             candMet = [candMet candMet(:,end)];
%                             candMet(:,6) = candMet(:,5);
%                             candMet(:,5) = focusMetric;
%                             candMet(:,4) = [];
%                             focusMetric((1-corrEllip)>0.3) = NaN;
%                             
%                             %Plane Consolidation occur here
%                             [part] = obj.planeConsolidation(candMet,focusMetric,fCandMet(:,3));
% 
%                             %we delete empty cells from the array
%                             idx2Empty = cellfun(@isempty,part);
%                             part(idx2Empty(:,1),:) = [];
%                    
%                             particleList{idx} = part;
%                             nParticles(idx) = length(part);
%                             
%                             if ~isempty(part)
%                                 idx2TP(idx) = idx;
%                             end
%                         
%                     end
%                     waitbar(i/nFrames,h,['Consolidating candidate... ' num2str(i) '/' num2str(nFrames) ' done']);
%                 end
%                 close(h);
%                 
%                 %#3 Storing List
%                 particle.List       = particleList;
%                 particle.nParticles = nParticles;
%                 particle.tPoint     = nFrames;
%                 particle.idx2TP     = nonzeros(idx2TP);
%                 particle.Traces     = [];
%                 particle.nTraces    = [];
%                 
%                 fileName = sprintf('%s%sparticle.mat',obj.raw.movInfo.Path,'\');
%                 save(fileName,'particle');
%             end
%             %#4 Storing particles in the object
%             obj.particles = particle;
%         end
        
        function showCorrLoc(obj)
            part = obj.particles.List;
            
            figure()
            hold on
            for i = 1:length(part)
                
                if ~isempty(length(part{i}))
                    
                    for j = 1: length(part{i})
                        data = part{i}{j}(3,1:3);
                        sizeMarker = 5;
                        scatter3(data(1),data(2),data(3),sizeMarker,data(3),'filled');
                       
                        
                    end
                end
                    
            end
            
            title('all Localization plotted');
            xlabel('x position in nm');
            ylabel('y position in nm');
            zlabel('z position in nm');
            
            
            hold off
            
        end
        
    end
    
    methods (Static)
        
        function [isPart]   = isPartPlane(current, next, fCand)
            %This function aim at determining whether a candidate from one
            %plane and the another are actually the same candidate on
            %different plane or different candidate. The decision is based
            %on threshold on localization distance, ellipticity and focus
            %metric.
            assert(size(current,2) == size(next,2), 'Dimension mismatch between the tail and partner to track');
            
            Thresh = 150; %~1px (100 nm) max as we SR-corrected the data 
            [checkRes1] = Core.MPParticleMovie.checkEuDist(current(:,1:2),...
                next(:,1:2),Thresh);
            
            % Test Z
            checkRes2 = fef;
            for i = 1: length(fCand)
                if fCand(i)
                threshold = 500; %nm
                [Res2] = Core.MPLocMovie.checkZPos(current(:,3),...
                next(i,3),threshold);
                else
                    Res2 = true;
                end
                checkRes2(i) = Res2;
            end
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
        
        function [res] = checkZPos(current,next,threshold)
            
            res = (abs(current-next)<threshold);
            
        end
        
    end
    
    
    methods (Access = private)
        
        function [corrData] = applyTrans(obj, data2Corr, transMat, refPlane, currentPlane)
            %act depending on whether the current plane is smaller or
            %bigger than the user-selected reference plane
            
            if currentPlane < refPlane

                idx2Corr = currentPlane:refPlane-1;
                sign = -1;

            elseif currentPlane > refPlane

                idx2Corr = refPlane:currentPlane-1;
                idx2Corr = fliplr(idx2Corr);
                sign = +1;

            else

                idx2Corr = [];

            end
            %1 Translation
            if ~isempty(idx2Corr)
                for j = 1:length(idx2Corr)
                    data2Corr = data2Corr + sign* transMat{idx2Corr(j)};
                end
            end
            
            corrData = data2Corr;
        
        end
        
        function [corrData] = applyRot(obj, data2Corr, rotMat, refPlane, currentPlane)
            
            %act depending on whether the current plane is smaller or
                %bigger than the user-selected reference plane
                if currentPlane < refPlane
                    
                    idx2Corr = currentPlane:refPlane-1;
                    sign = false;
                    
                elseif currentPlane > refPlane
                    
                    idx2Corr = refPlane:currentPlane-1;
                    idx2Corr = fliplr(idx2Corr);
                    sign = true;
                    
                else
                    
                    idx2Corr = [];
                    
                end
                
                %Rotation
                if ~isempty(idx2Corr)
                    %Pad Z coordinate
                    data2Corr(:,3) = 0;
                    %remove center of mass (CM)
                    CM = mean(data2Corr);
                    data2Corr = data2Corr - CM;
                    
                    %Change the orientation of the data (should be [x;y;z]
                    %not [x y z]
                    data2Corr =  data2Corr';
                    
                    %Correction occur here
                    for j = 1:length(idx2Corr)
                        
                        if sign
                            rot = rotMat{idx2Corr(j)}';
                        else
                            rot = rotMat{idx2Corr(j)};
                        end
                           
                            data2Corr  = (rot*data2Corr);
                            data2Store = data2Corr';
                            corrData   = data2Store(:,1:2)+CM(1:2);
                    end
                else
                    corrData = data2Corr(:,1:2);
            
                end
  
        end
    
        function [zPos,inRange] = getZPosition(obj,ellip,zCal,currentPlane)
            
            relZ = obj.calibrated.oRelZPos;
                       
            zVec = -2000:1:2000; %Here we assume accuracy >= 1nm
                       
            fit = polyval(zCal.calib{currentPlane,1},zVec);
            %fit = ppval(zCal{particle(3,end),2},zVec);%spline
            
            %find the index of the value the closest to the particle
            %ellipticity
            
             [~,idx] = min(abs(fit-ellip));
             
             zPos = zVec(idx)+ relZ(currentPlane)*1000;          
             inRange = and(ellip>=zCal.fitZParam(1).ellipRange(1),...
                 ellip<= zCal.fitZParam(1).ellipRange(2));

        end
        
        function [zRange] = getZRange(obj,ellipRange,zCal,currentPlane)
            
            relZ = obj.calibrated.oRelZPos;
                       
            zVec = -2000:1:2000; %Here we assume accuracy >= 1nm
                       
            fit = polyval(zCal.calib{currentPlane,1},zVec);
            %fit = ppval(zCal{particle(3,end),2},zVec);%spline
            
            %find the index of the value the closest to the particle
            %ellipticity
            
             [~,idx1] = min(abs(fit-ellipRange(1)));
             [~,idx2] = min(abs(fit-ellipRange(2)));
             
             zPos1 = zVec(idx1)+ relZ(currentPlane)*1000;    
             zPos2 = zVec(idx2)+ relZ(currentPlane)*1000;
             
             zRange = [zPos1, zPos2];
             
        end
        
        function [candidateList] = planeConsolidation(obj,candMet,focusMetric,fCandMet)
            %Loop through all candidate of a given frame and match them
            %between frame until none can be match or all are matched.
            nPlanes = obj.calibrated.nPlanes;
            counter = 1;
            nPart = 0;
            maxIt = size(candMet,1);
            candidateList = cell(max(size(find(~isnan(focusMetric)))),1);
            ellipRange = obj.ZCal.fitZParam.ellipRange;
            %continue until the list is not empty
            while and(~isempty(focusMetric), ~isnan(nanmax(focusMetric)))
                
                if counter> maxIt
                    
                    error('While loop ran for an unexpectedly long time, something might be wrong');
                    
                end
                
                %Find candidate in best focus
                [~,idx] = max(focusMetric);
                currentCand = candMet(idx,:);
                currentPlane = candMet(idx,end);
                
                %Check which planes are to be checked (currently 2 planes
                %above and 2 planes below the given plane
                planes2Check = currentPlane-2:currentPlane-1;
                planes2Check = planes2Check(planes2Check>0);
                planes2Check = [planes2Check currentPlane+1:currentPlane+2];
                planes2Check = planes2Check(planes2Check<nPlanes+1);
                
                particle = nan(5,6);
                particle(3,:) = currentCand;
                nCheck = length(planes2Check);
                for i = 1:nCheck
                    
                    cand = candMet(candMet(:,end) == planes2Check(i),:);
                    fCand = fCandMet(candMet(:,end) == planes2Check(i),:);
                    [isPart] = Core.MPLocMovie.isPartPlane(currentCand,cand,fCand);
                    if ~all(isPart ==0)
                        id = cand(isPart,end)-currentCand(end);
                        particle(round(length(currentCand)/2)+id,:) = cand(isPart,:);
                    end
                    
                end
                
                %Check if the resulting configuration of the plane make
                %sense e.g. no hole in the configuration
                planeConfig = particle(:,end);
                [checkRes] = Core.MPParticleMovie.checkPlaneConfig(planeConfig,nPlanes);
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