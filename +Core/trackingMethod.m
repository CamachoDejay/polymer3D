classdef trackingMethod < handle
    %"Static class" to hold method related to tracking
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
 
    end
    
    methods(Static)
        function listIdx  = connectParticles(List,listBool,idx, trackParam)
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
                nPartInFrame = length(checkRes);
                
                if ~all(checkRes==0)
                    %We use reshape to input the different particles of the next
                    %frame at once by storing them in the 3rd dimension
                    nextPart = List{currentIdx(1)+1};
                    nextPart = nextPart(logical(checkRes));
                    
                    shape = [size(nextPart{1},1) size(nextPart{1},2) nPartInFrame];%Unused particle are expected to be still 5x5
                    nextParts(:,:,1:nPartInFrame) = reshape(cell2mat(List{currentIdx(1)+1}'),shape);
                    nextParts = nextParts(:,:,logical(checkRes));
                    [isPart] = Core.trackingMethod.isPartFrame(part2Track,nextParts,1,trackParam);
                    
                   
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
                roughcheck1 = squeeze(abs(current(3,end)-next(3,end,:))<=1);
                
                if all(roughcheck1 ==0)
                    disp('Something is wrong your focus changed more than one plane between 2 frames');
                else
                    %Check that at least 2 planes are in common
                    commonPlanes = Core.trackingMethod.findCommonPlanes(current(:,end),squeeze(next(:,end,:)));
                    roughcheck2  = sum(commonPlanes,1)>=nConsistentPlanes;
                    if all(roughcheck2 ==0)
                        disp('Less than 2 planes in common, breaking out');
                    else
                        
                        for i = 1 : size(next,3)
                            % Test Euclidian distance
                            Thresh = trackParam.euDistPx; %in px
                            [checkRes1] = Core.MPParticleMovie.checkEuDist(current(commonPlanes(:,1,i),1:2),...
                                squeeze(next(commonPlanes(:,2,i),1:2,i)),Thresh);
                            
                            % Test ellipticity
                            eWeight = [1 2 3 2 1];
                            thresh = trackParam.ellip;
                            [checkRes2] = Core.MPParticleMovie.checkEllipticity(current(commonPlanes(:,1,i),3),...
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
end

