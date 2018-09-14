classdef MPTrackingMovie < Core.MPLocMovie
    %trackingMovie will take particles detected and try to match them along
    %frames
       
    properties
        traces3D
    end
    
    methods
        function obj = MPTrackingMovie(raw, MPCal, SRCal, zCal)
            %trackingMovie Construct an instance of this class
            %   Detailed explanation goes here
             obj  = obj@Core.MPLocMovie(raw,MPCal,SRCal,zCal);
             
        end
        
        function trackParticle(obj,trackParam)
            
             %track the particle in the Z direction (3rd dimension here)
            assert(~isempty(obj.calibrated),'Data should be calibrated to do ZzCalibrationration');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before zzCalibrationration');
            assert(~isempty(obj.particles), 'No particles found, please run superResConsolidate method before doing ZzCalibrationration');
            assert(isstruct(trackParam),'Tracking parameter is expected to be a struct with two field euDistXY and euDistZ')
            assert(and(isfield(trackParam,'euDistXY'),isfield(trackParam,'euDistZ')),...
                'Tracking parameter is expected to be a struct with two field "euDistPXY" and "euDistZ"')
            assert(~isempty(obj.corrected),'Data needs to be corrected before tracking');
            assert(and(obj.corrected.XY,obj.corrected.Z),'Data needs to be corrected before tracking');
            
            %We copy the List as boolean to keep track of where there are
            %still particles left
            [listBool] = Core.trackingMethod.copyList(obj.particles.List,1);
            %We copy as NaN for storage of the traces;
            [traces]   = Core.trackingMethod.copyList(obj.particles.List,NaN);
            %We pick the first particle available
            [idx] = Core.trackingMethod.pickParticle(listBool);
            counter = 1;
            errCount =1;
            while (idx)
                %loop until there is no particle (pickParticle return false)
                if errCount>1000
                    warning('While loop ran for unexpectedly longer time');
                    break;
                    
                end
                %Connect particles (cf consolidation but across frames
                [listIdx] = Core.MPTrackingMovie.connectParticles(obj.particles.SRList,listBool,idx, trackParam);
                %if the particle was connected in less than 5 frames we remove
                % its appearance from the list bool
                if length(listIdx) < 5
                    
                    [listBool] = Core.trackingMethod.removeParticles(listBool,listIdx);
                    
                else
                    %Otherwise we store traces, increment counter and remove.
                    [traces]  = Core.trackingMethod.storeTraces(traces,listIdx,counter);
                    counter = counter +1;
                    [listBool] = Core.trackingMethod.removeParticles(listBool,listIdx);
                    
                end
                % We pick a new particle and start all over again
                [idx] = Core.trackingMethod.pickParticle(listBool);
                errCount = errCount +1;
            end
            counter = counter -1;
            
            obj.particles.traces = traces;
            obj.particles.nTraces = counter;
            
            [trace3D] = obj.get3DTraces;
            
            obj.traces3D = trace3D;
            
        end
        
        function showTraces(obj)
            traces = obj.traces3D;
            obj.showCorrLoc;
            
            gcf;
            
            hold on
            
            for i = 1: length(traces)
                currentTrace = traces{i};
                data = table2array(currentTrace(:,{'row','col','z'}));
                plot3(data(:,1), data(:,2), data(:,3))
                
            end
            hold off
            
            figure
            hold on 
            for i = 1: length(traces)
                currentTrace = traces{i};
                plot3(currentTrace.row, currentTrace.col, currentTrace.z)
                
            end
            hold off
        end
         
    end
    
    methods (Static)
        
        function [listIdx] = connectParticles(List,listBool,idx,trackParam)
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
                [checkRes] = Core.trackingMethod.checkListBool(listBool,currentIdx(1)+1);
                
                
                if ~all(checkRes==0)
                    %We use reshape to input the different particles of the next
                    %frame at once by storing them in the 3rd dimension
                    nextPart = List{currentIdx(1)+1};
                    nextPart = nextPart(logical(checkRes));
                    %used to be additional use of checkRes, why?
                    [isPart] = Core.MPTrackingMovie.isPartFrame(part2Track,nextPart,trackParam);
                    
                   
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
                        counter = counter+1;
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
        
        function [isPart]   = isPartFrame(current, next, trackParam)
            %This function is designed to have PSFE plate ON
            assert(and(istable(current),iscell(next)), 'unexpected format in partners to track');
            
                %ZStack, consolidation between frame
                %The calculation here is ran in parallel, we check if
                %the current particle is a partner of one of the
                %particles in the next frame. Some indexing or step
                %might therefore seems unnecessary but allow to take
                %any number of particles
   
                nPart = length(next);
                isPart = zeros(nPart,1);
                
                for i = 1 : nPart
                    
                    nextPart = next{i};
                    
                    % Test Euclidian distance
                    ThreshXY = trackParam.euDistXY; %in nm
                    ThreshZ  = trackParam.euDistZ;
                    [checkRes1] = Core.MPParticleMovie.checkEuDist([current.row current.col],...
                        [nextPart.row, nextPart.col],ThreshXY);
                    
                    [checkRes2] = Core.MPTrackingMovie.checkZ(current.z,nextPart.z,ThreshZ);
                    
                    %Both test need to pass to be partenaires
                    isPart(i) = checkRes1.*checkRes2;
                           
                end
%                 test = and(length(isPart)>1,all(isPart==1));
%                 
%                 %dbstop at 173 if test==1
%                 
                isPart = logical(isPart);
                
            
        end
        
        function [checkRes] = checkZ(Z1,Z2,Thresh)
            
            checkRes = abs(Z1-Z2) <= Thresh;
        end 
        
                 
    end
    
    methods (Access = private)
        
        function [traces3D ] = get3DTraces(obj)
            partList = obj.particles.SRList;
            traces = obj.particles.traces;
            nTraces = obj.particles.nTraces;
            
            traces3D = cell(nTraces,1);
            fCounter = ones(nTraces,1);
            for i = 1 : length(partList)
                
                currentParts = partList{i};
                currentTraces = traces{i};
                for j = 1 : length(currentParts)
                    
                    currentPart = currentParts{j};
                    currentTrace = currentTraces{j};
                    
                    if ~isnan(currentTrace)
                        
                        traces3D{currentTrace}(fCounter(currentTrace),...
                            {'row','col','z','frame'}) = [currentPart(:,{'row','col','z'}),  table(i)];
                        
                        
                        fCounter(currentTrace) = fCounter(currentTrace)+1;
                    else
                    end
                    
                
                end
            end
            
            
            
        end
        
       
        
    end
end

