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
                
                
                if ~all(checkRes==0)
                    %We use reshape to input the different particles of the next
                    %frame at once by storing them in the 3rd dimension
                    nextPart = List{currentIdx(1)+1};
                    nextPart = nextPart(logical(checkRes));
                    %used to be additional use of checkRes, why?
                    [isPart] = Core.trackingMethod.isPartFrame(part2Track,nextPart,1,trackParam);
                    
                   
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
            commonPlanes(:,1) = ismember(planeInCurrent,planeInNext);
            commonPlanes(:,2) = ismember(planeInNext,planeInCurrent);

            commonPlanes = logical(commonPlanes);
        end
        
        function [isPart]   = isPartFrame(current, next, direction, trackParam)
            %This function is designed to have PSFE plate ON
            assert(abs(direction) == 1, 'direction is supposed to be either 1 (up) or -1 (down)');
            assert(and(istable(current),iscell(next)), 'unexpected format in partners to track');
            
                %ZStack, consolidation between frame
                %The calculation here is ran in parallel, we check if
                %the current particle is a partner of one of the
                %particles in the next frame. Some indexing or step
                %might therefore seems unnecessary but allow to take
                %any number of particles
                
                nConsistentPlanes = trackParam.commonPlanes;
                nPart = length(next);
                isPart = zeros(nPart,1);
                
                for i = 1 : nPart
                    
                    nextPart = next{i};
                    %check focus is not more than one plane away
                    roughcheck1 = squeeze(abs(current.plane(3)-nextPart.plane(3))<=2);

                    if roughcheck1 ==0
                        %disp('Something is wrong your focus changed more than one plane between 2 frames');
                    else
                        %Check that at least 2 planes are in common
                        commonPlanes = Core.trackingMethod.findCommonPlanes(current.plane,nextPart.plane);
                        roughcheck2  = sum(commonPlanes,1)>=nConsistentPlanes;
                        if all(roughcheck2 ==0)
                            disp('Less than 2 planes in common, breaking out');
                        else

                                % Test Euclidian distance
                                Thresh = trackParam.euDistPx; %in px
                                [checkRes1] = Core.MPParticleMovie.checkEuDist([current.row(commonPlanes(:,1)) current.col(commonPlanes(:,1))],...
                                    [nextPart.row(commonPlanes(:,2)), nextPart.col(commonPlanes(:,2))],Thresh);

                                % Test ellipticity
                                [checkRes2] = Core.MPParticleMovie.checkEllipticity(current.ellip(commonPlanes(:,1)),...
                                    nextPart.ellip(commonPlanes(:,2)),direction);
                                
                                %To be a particle, we want the position and ellipticity to be
                                %consistent in at least 2 planes 
                                isPart(i) = and(length(find(checkRes1))>=nConsistentPlanes,...
                                    length(find(checkRes2))>=nConsistentPlanes);

                        end
                    end
                end
                
                isPart = logical(isPart);
                
            
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
        
        %%%%%%%%%%%%%%%%%%% FROM SERGEY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [ALLDataConverted,AllFieldName] =ConvertData(input_data,ImMax)
    
            assert(GetInitialDataType(input_data)==0,'Wrong format of the input data :(. Can not continue, sorry for the inconvenience!')
            AllFieldName = fieldnames(input_data);
            AllFieldName(end-2:end) = [];
            AllFieldName{end+1} = 'NumParticle';
            input_data = table2array(input_data);
            ALLDataConverted = cell(ImMax,1);

            for i=1:ImMax         
                    ALLDataConverted{i} =input_data(input_data(:,end)==i,:);
            end

        end

        
        %%% RESOLVE TRACKING CONFLICTS
        function AllCandidates = ResolveConflicts(InitialArray)
              % Preallocate
              %allarray = zeros( size(InitialArray,2),3);
              %Boris edit: Correct preallocation
              totSize = cell2mat(cellfun(@size,InitialArray(2,:),'UniformOutput',false));
              totSize = sum(totSize(1:2:end));
              allarray = zeros(totSize,3);
              
              % Reform array into Column 1: Candidate indeces in NextFrame
              %                    Column 2: Distance to the point particle id index in
              %                    previous frame
              %                    Column 3: Corresponding particle-id indeces in previous frame
              for i=1:length(InitialArray)
                 InitialArray{2,i}(:,3)=i;
                 [m,~]=size(InitialArray{2,i});
                 ind = find(allarray(:,1)==0,1,'first');
                 allarray(ind:ind+m-1,:) = [InitialArray{2,i}];
              end

              % Find the overall minima of the array sum using Hungarian algorithm, and
              % add it to the AllCandidates array of all of possible further
              % tracked particle trajectories 
              allarray(allarray(:,1)==0,:)=[];      

              dim = max([max(allarray(:,1)),max(allarray(:,3))]);
              MatrixToMinimize = Inf(dim,dim);
              for i=1:length(allarray(:,1))
                MatrixToMinimize(allarray(i,1),allarray(i,3))=allarray(i,2);
              end
              [AllCandidates,~] = Minimization.munkres(MatrixToMinimize) ;
              AllCandidates(2,:)=1:length(AllCandidates);
              AllCandidates(:,AllCandidates(1,:)==0) =[];
              AllCandidates = fliplr(AllCandidates');


        end

        %%% FIND POSSIBLE NEIGHBOURS
        function PossibleNeighbours = FindPossibleNeighbours(NextFrame,PreviousFrame,Radius,PossibleNeighbours)
               % Searches for any possible neighbours in the radius = Radius. Outputs a cell array of indeces from NextFrame, corresponding to any points 
               %  withing the circle of radius = Radius. Removes the used line from
               %  the array every iteration step, and stops once array is empty
               if ~isempty(PreviousFrame)

                   DistanceToAllThePointsInTheNextFrame = Core.trackingMethod.dist3d(PreviousFrame(1,1:3), NextFrame(:,1:3));
                   [DistancesSmallerThenRadius] = find( DistanceToAllThePointsInTheNextFrame < Radius );

                   PossibleNeighbours(end+1)  = {[DistancesSmallerThenRadius,DistanceToAllThePointsInTheNextFrame(DistancesSmallerThenRadius)]};
                   PossibleNeighbours = Core.trackingMethod.FindPossibleNeighbours(NextFrame,PreviousFrame(2:end,:),Radius,PossibleNeighbours);  

               end
               NoParticlesDetected = 0; 

               for i=1:length(PossibleNeighbours)
                   if ~isempty(PossibleNeighbours{1,i})
                      NoParticlesDetected = 1; 
                   end
               end
               if  logical(NoParticlesDetected)==0
                    warning('No neighbours detected in the next frame, check your radius or detection threshold')
               end
        end

        %%% FIND NEIGHBOURS IN NEXT FRAME
        function PossibleNeighbours = SearchNeighbours(input_data ,IMPREV,radius)            
            allids = input_data.ConstructACellArrayWithAllParticleIds(IMPREV,{[]});
            tempneighbours = Core.trackingMethod.FindPossibleNeighbours(input_data.dataNext,IMPREV,radius,{[]});
            PossibleNeighbours =[allids(2:end) ; tempneighbours(2:end)];
        end
        %%% 3D EUCLIDIAN DISTANCES
        function distan = dist3d(Point, Arr)

            distan = sqrt((Arr(:,1)-Point(1)).^2 + (Arr(:,2)-Point(2)).^2 + (Arr(:,3)-Point(3)).^2);
        end
        
        %%% This function adds new tracked data to the already existing array. If new particles come in, it will expand the array and maxid as well.
        function [DataAdded,maxid] = AddDataToTracked(PreviousData,DataToAdd, NeighbourData,maxid)
        %Preallocate and initialize
           DataAdded =[];

           if ~isempty(NeighbourData)
           [~,n]= size([DataToAdd(NeighbourData(1,1),:),PreviousData(NeighbourData(1,2),end)]);
           DataAdded =zeros(1000,n);

        %Add tracked data

           if ~isempty(DataToAdd)
               CopyOfNeighbourData= NeighbourData;
               DataAdded= [DataToAdd(NeighbourData(:,1),:) ,PreviousData(NeighbourData(:,2),end)];
        %Remove all already assigned ids from the datatoadd, so that only unassigned ids remain          
               DataToAdd(CopyOfNeighbourData(:,1),:)=[];
        %Run through unassigned data and add as new id, and increase maxid by 1
               if ~isempty(DataAdded)
               newmaxid =maxid+length(DataToAdd(:,1));
               DataAdded =[DataAdded ; [DataToAdd ,[maxid+1:newmaxid]']];
               maxid = newmaxid;
               end

           else
           warning('No Particles were detected , check your radius threshold')
           return
           end
           DataAdded(DataAdded(:,1)==0,:)=[];
           end
        end
%% FUNCTION DEFINITIONS FOR MEMORY HANDLING
%%% Cleans memory from irrelevant points for which time has exceed the maximum time allowed for memory.
        function CleanedMemory = CleanMemory(data, Time,MaxTime)
            if ~isempty(data)   

               data(abs(data(:,end-1)-Time)>MaxTime,:) = [];
               CleanedMemory = data;
            else
               CleanedMemory = [];
            end


        end
        %%% Adds to memory any particles that did not find a neighbour in their vicinity
        function AddedMemory = AddToMemory(Neighbours, PreviousData)
            PreviousData(Neighbours(:,2),:) = [];
            AddedMemory = PreviousData;
        end


%% Conversion of final output

        function Data = ConvertFinalOutput( TrackedData ,Data,AllFieldName)
    
            while ~isempty(TrackedData)
                if ~isempty(TrackedData{1})
                TimeFrame =TrackedData(1);
                IdMax = max(TimeFrame{1,1}(:,end));
                for i=1:IdMax
                    Data{i} = [Data{i} ;TimeFrame{1,1}(TimeFrame{1,1}(:,end)==i,:)];
                end
                end
                TrackedData(1) =[];
            end
            for i=1:length(Data)

               Data{i} = array2table(Data{i}); 
               Data{i}.Properties.VariableNames = AllFieldName;
            end

        end        
    end
end

