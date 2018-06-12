function particleMetric = posConsolidation(data,candidate)
delta = 6;
candidate = {candidate};
sig = [108.33/95 108.33/95];
particleMetric = cell(length(candidate),1);
for i = 1 : length(candidate)
    
    frameCandidate = candidate{i};
    currentk = 1;
    %probably need to call getFrame somewherehere and getCandidate
    partMet = zeros(length(frameCandidate),7);
    for j = 1 : numel(fields(data))
        planeCandidate = frameCandidate(frameCandidate(:,3)==j,1:2);
        planeData = data.(sprintf('plane%d',j));
        
        pM = zeros(length(planeCandidate),6);
        for k = 1 : length(planeCandidate)
        [roi_lims] = EmitterSim.getROI(planeCandidate(k,2), planeCandidate(k,1),...
            delta, size(planeData,2), size(planeData,1));
        ROI = planeData(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2));
        
        [x,y,e] = Localization.phasor(ROI);
        [LRT,RMSD] = Localization.likelihoodRatioTest(ROI,sig,[y x]);
        xPos = planeCandidate(k,2) + x;%Check with Rafa
        yPos = planeCandidate(k,1) +y;
        %COMMENT COLUMN MEANING
        pM(k,:) = [xPos yPos e LRT RMSD k];
        end
        
        partMet(currentk:currentk+k-1,:) = [pM j*ones(k,1)];
        currentk = currentk+k;
    end
    counter = 1;
    nPart = 0;
    maxIt = size(partMet,1);
    while ~isempty(partMet)
        if counter> maxIt
            error('While loop ran for an unexpectedly long time, something might be wrong');
        end
        
        %Find candidate in best focus
        [~,idx] = max(abs(partMet(:,4)));
        
        if abs(partMet(idx,3)-1) > 0.2 %Discuss with Rafa 
            idxF = idx;
        else
            nPart = nPart +1;
            %Check if there are planes up
            [idx1] = checkPlane(partMet,idx,1);
            
            %Check if they are planes down
            [idx2] = checkPlane(partMet,idx,-1);
            
            idxF = [idx1 idx idx2];
            
            %Store
            particleFinalList.(sprintf('particle%d',nPart)) = partMet(idxF(:),:);
            %Remove from list
            
        end
        %We remove the particle(s) from the list
        partMet(idxF,:) = [];
        
        counter = counter+1;
    end
    
    particleMetric{i} = particleFinalList;
    
    
end
end

function [idx2Part] = checkPlane(partMet,idx,direction)
 %This function is designed to work on 8 imaging planes only.
 nP = 8; %TODO: extract from object
assert(abs(direction) == 1, 'direction is supposed to be either 1 (up) or -1 (down)');
%assert partMet
currentPlane = partMet(idx,7);
currentPart  = partMet(idx,:);
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
    partInPlane = partMet(partMet(:,end) == planeIdx,:);
    %Check position %Discuss with Rafa
    partInPlane = partInPlane(and(abs(partInPlane(:,1)-currentPart(1))<2,...
                                  abs(partInPlane(:,2)-currentPart(2)<2)),:);
    if size(partInPlane,1)>1
        disp('More than one particle in close proximity')
    end
    if(isempty(partInPlane))
        disp('No other candidate in plane with similar position')
    else
    %Check ellipticity (only if PSFE was on
        switch direction
            case 1
                partInPlane = partInPlane(partInPlane(:,3) > currentPart(3),:);
            case -1
                partInPlane = partInPlane(partInPlane(:,3) < currentPart(3),:);
        end
    
        if(isempty(partInPlane))
            disp('No close candidate in plane consistent with ellipticity')
        else
            %check GLRT
            partInPlane = partInPlane(abs(partInPlane(:,4)) < abs(currentPart(4)),:);
            if(isempty(partInPlane))
            disp('No close candidate in plane consistent with the LRT parameter')
            end
        end
    end
    
    if ~isempty(partInPlane)
    newIdx = find(partMet(:,1) == partInPlane(1));
    idx2Part = [idx2Part newIdx];
    end
end

end


