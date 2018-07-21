classdef trackingMovie < Core.mpLocMovie
    %trackingMovie Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        traces
        zCalibration
    end
    
    methods
        function obj = trackingMovie(raw, cal, zCalibration)
            %trackingMovie Construct an instance of this class
            %   Detailed explanation goes here
            obj  = obj@Core.mpLocMovie(raw,cal);
            
            obj.zCalibration = zCalibration;
        end
        
        function set.zCalibration(obj, zCalibration)
            assert(isstruct(zCalibration),'zCalibration is expected to be a struct with 4 fields')
            assert(and(and(isfield(zCalibration,'file'),isfield(zCalibration,'data')),...
                and(isfield(zCalibration,'fitZParam'),isfield(zCalibration,'trackParam'))),...
                'zCalibration is expected to be a struct with 4 fields');
            
            obj.zCalibration = zCalibration;
            
        end
        
        function giveZCal(obj,zCalibration)
            
            obj.zCalibration = zCalibration;
            
        end
        
        function [traces] = get3DTraces(obj)
            %Extract 3D traces for each particles
            list = obj.particles.List;
            tracesIdx = obj.particles.Traces;
            pxSize = obj.info.pxSize;
            %ellipt range used for fitting
            elliptRange = [obj.zCalibration.syncEllip{1,3}(1) obj.zCalibration.syncEllip{1,3}(2)];
            elliptRange = elliptRange(1):0.01:elliptRange(2);
            %we weigh the average later base on how much out of focus the
            %plane was.
            wRange1 = length(elliptRange(elliptRange<=1));
            wRange2 = length(elliptRange(elliptRange>=1));
            weight1 = linspace(1,5,wRange1);
            weight2 = linspace(5,1,wRange2);
            finalWeight = [weight1 weight2];
            
            if isempty(obj.zCalibration.cal)
                
                warning('no z zCalibrationration detected, only show 2D plot');
                
            end
            %Calculation of x-y-z position occurs within the loop here
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
                        %Weighed average
                        xAvg = sum(diag(list{i}{j}(id,2)* weight))/sum(weight) * pxSize;
                        yAvg = sum(diag(list{i}{j}(id,1)* weight))/sum(weight) * pxSize;
                        %best focus Value
                        x = list{i}{j}(3,2)* pxSize;
                        y = list{i}{j}(3,1)* pxSize;
                        
                        if ~isempty(obj.zCalibration.cal)
                            %Calculation for Z is occur here
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
            %Storage
            obj.particles.traces3D = traces;
        end
        
        function [zPos,zAvg] = getZPosition(obj,particle,elliptRange,finalWeight)
            %extract Z position based on ellipticity and zCalibrationration curve
            assert(~isempty(obj.zCalibration),'No zCalibration found, please give the zCalibration using giveZCal before calculating Z position');
            zCal = obj.zCalibration.cal;
            relZ = obj.calibrated.oRelZPos;
            syncEllip = obj.zCalibration.data;

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
            %Calculate z occur withing the loop
            for k = 1 :size(data2Keep,1)

            [~,idx1(k)] = min(abs(elliptRange-data2Keep(k,3)));

            zVectmp = min(syncEllip{data2Keep(k,end),1}) : 1 : max(syncEllip{data2Keep(k,end),1});
            fit = polyval(zCal{data2Keep(k,end),1},zVectmp);

            [~,idx] = min(abs(fit-data2Keep(k,3)));
            zAvg(k) = zVectmp(idx) + relZ(data2Keep(k,end))*1000;
            end

            weight = finalWeight(idx1);
            %weighed average occur here
            %diag is taken because of the matrix operation between weight
            %and zAvg
            zAvg = sum(diag(zAvg(:)* weight(:)'))/sum(weight);

        end
        
    end
end

