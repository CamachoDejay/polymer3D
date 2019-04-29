classdef MPZCalMovie < Core.MPCalMovie
    %zCal will hold the information related to zCalibration movies as well as all
    %the method linked to the zzCalibrationration.
    
    properties (SetAccess = 'protected')
        zData
       
    end
    
    methods
        
        function obj = MPZCalMovie(raw,cal,info)
            %UNTITLED6 Construct an instance of this class
            %   Detailed explanation goes here
            obj  = obj@Core.MPCalMovie(raw,cal,info);
            
        end
        
        function set.zData(obj, zData)
            
            obj.zData = zData;
            
        end
        
        function [zCalData] = getCalData(obj,traces,nPart)
            %Extract all the data across frame separated by planes (1zCalibration
            %curve by plane
            zCalData = cell(obj.calibrated.nPlanes,nPart);
            
            for i = 1:length(traces)
                
                if isempty(traces{i})
                    
                else
                    
                    for j = 1:length(traces{i})
                        if ~isnan(traces{i}{j})
                            currentParticle = obj.particles.List{i}{j}; 
                            planes = currentParticle.plane;
                            planes = planes (~isnan(planes));
                            
                            for k = 1 : length(planes)
                                data2Store = obj.particles.List{i}{j};
                                data2Store = data2Store(data2Store.plane == planes(k),:);
                                data2Store.frame = i;
                                zCalData{planes(k),traces{i}{j}}(i,:) = data2Store ;
                                
                            end
                        else
                        end
                        
                    end
                end
            end
            %clean up the data
            for i = 1 : size(zCalData,1)
                for j = 1 : size(zCalData,2)
                    
                    if ~isempty(zCalData{i,j})
                        
                        zCalData{i,j}(zCalData{i,j}.row==0,:) = [];
                        
                    else
                    end
                    
                end
            end
            obj.zData.calData = zCalData;
        end
        
        function showParticlesTracked(obj,ips)
            %display a movie where both localization and consolidation are
            %showed. The display needs to be improved !
            %TODO : improve display of the movie
            assert(~isempty(obj.traces),'You need to get the traces before displaying them');
            
            nPlanes = obj.calibrated.nPlanes;
            colors = rand(obj.traces.nTrace,3);
            
            for i = 1 : length(obj.traces.trace)
                nParticles = length(obj.traces.trace{i});
                obj.showCandidate(i);
                
                h = gcf;
                
                for j = 1 : nPlanes
                    subplot(2,nPlanes/2,j)
                    hold on
                    for k = 1 : nParticles
                        currPart = obj.particles.List{i}{k};
                        labelColor = obj.traces.trace{i}{k};
                        if ~isnan(labelColor)
                            if(~isempty(currPart(currPart(:,end) == j)))
                                part2Plot = currPart(currPart(:,end) == j,:);
                                plot(part2Plot(2),part2Plot(1),'o',...
                                    'LineWidth',2, 'MarkerSize',10, 'MarkerEdgeColor',colors(labelColor,:));
                                drawnow
                            end
                        end
                    end
                    hold off
                    pause(1/ips);
                end
            end
        end
        
        function [zSyncCalData] = syncZCalData(obj,zCalData,fitZParam)
            %Fit the ellipiticty - zPos data for each particles and
            %synchronized them by putting the point where ellipticity = 1
            %to zPos=0+PlanePosition. It uses the zCaldata which represent
            %the data point for each particle for a particle plane. Fitting
            %is done by polynomial to increase accuracy (if e = 1.05 then
            %this point will be locate not exactly at 0 but will be shifted
            assert(isstruct(fitZParam), 'fitZParam is expected to be a struct with two fields, deg for polynomial fit and ellipRange for ellipticity to consider');
            assert(and(isfield(fitZParam,'deg'),isfield(fitZParam,'ellipRangeCal')),...
                'fitZParam is expected to be a struct with two fields, deg for polynomial fit and ellipRange for ellipticity to consider');
            
            nPlanes = obj.calibrated.nPlanes;
            deg = fitZParam.deg;
            minEllipt = fitZParam.ellipRangeCal(1);
            maxEllipt = fitZParam.ellipRangeCal(2);
            [zStep,~] = obj.getZPosMotor;
            zStep = zStep(1);
            zSyncCalData = cell(nPlanes,2);
            
            for i = 1:size(zCalData,1)
                for j = 1:size(zCalData,2)
                    %Extract data in acceptable ellipticity range
                    currentData = zCalData{i,j};
                    if ~isempty(currentData)
                        frames = currentData.frame;
                        zPos = frames*zStep;
                        fMetric = currentData.fMetric;
                        ellipt = currentData.ellip;
                        zPos = zPos(and(ellipt<maxEllipt, ellipt > minEllipt));
                        fMetric = fMetric(and(ellipt<maxEllipt, ellipt > minEllipt));
                        frames = frames(and(ellipt<maxEllipt, ellipt > minEllipt));
                        ellipt = ellipt(and(ellipt<maxEllipt, ellipt > minEllipt));
                        
                         %Let us filter bad data out:
                        lfMet = mean(fMetric)-std(fMetric);
                      
                        fMetric = fMetric(fMetric>lfMet);
                        frames = frames(fMetric>lfMet);
                        ellipt = ellipt(fMetric>lfMet);
                        zPos = zPos(fMetric>lfMet);
                        
                        %Now we shift in z the value closest to ellipt =1 for
                        %the fit == rough synchronization
                        [~,idx] = min(abs(ellipt-1));
                        focus1 = zPos(idx);
                        zPos = zPos - focus1;
                        
                        if length(ellipt)<=deg
                            %If not enough data for good fit we do not do
                            %anyhting
                            %we only process the data if there are point above
                            %and below focus.
                        elseif and(~isempty(ellipt(ellipt<1)),~isempty(ellipt(ellipt>1)))
                            
                            p = spline(zPos,ellipt);
                            zVec = min(zPos):0.001:max(zPos);%1nm step
                            fit = ppval(p,zVec);
                            %p = polyfit(zPos,ellipt,deg);
                            
                            %Fitting occurs here
                           % fit = polyval(p,zVec);
                            
    
                            [~,idx] = min(abs(fit-1));
                            focus2 = zVec(idx);
                            shift = focus1+focus2;%take into account both synchronization
                            
                            zToNm = (currentData.frame*zStep - shift)*1000;
                            ellip2Store = currentData.ellip;
                            magX2Store = currentData.magX;
                            magY2Store = currentData.magY;
                            data2Store = table(zToNm,ellip2Store,'VariableNames',{'z','ellip'});
                            data2Store.magX = magX2Store;
                            data2Store.magY = magY2Store;
                            data2Store.fMetric = zCalData{i,j}.fMetric;
                            data2Store.gFitMet =  zCalData{i,j}.gFitMet;

                            %tmp Store
                            zSyncCalData{i,1} = [zSyncCalData{i,1}; data2Store];
                            zSyncCalData{1,2} = [zSyncCalData{1,2}; data2Store];

                        end
                    end
                end
             end
            %tmp Store
            if ~isempty (zSyncCalData{i,1})
                [~,ind] = sort(zSyncCalData{i,1}.z);
                zSyncCalData{i,1} = zSyncCalData{i,1}(ind,:);
            end
                
            
            %final storing and output
            if ~isempty(zSyncCalData{1,2})
                [~,ind] = sort(zSyncCalData{1,2}.z);
                zSyncCalData{1,2} = zSyncCalData{1,2}(ind,:);
            end
            zSyncCalData{1,3} = [minEllipt, maxEllipt, deg];
            obj.zData.syncEllip = zSyncCalData;
        end
        
        function [traces,counter] = zTracking(obj,trackParam)
            %track the particle in the Z direction (3rd dimension here)
            %Here we do not expect any big movement from one frame to the
            %other so we give a warning if the tracking parameter seems to
            %soft.
            assert(isstruct(trackParam),'Tracking parameter is expected to be a struct with two field "euDistPx" and "ellip"')
            assert(and(isfield(trackParam,'euDistPx'),isfield(trackParam,'ellip')),...
                'Tracking parameter is expected to be a struct with two field "euDistPx" and "ellip"')
            
            if trackParam.euDistPx > 1
                warning('Current euclidian distance thresholds is probably too high, we do not expect much movement from one frame to the next here')
            end
            
            if or(trackParam.ellip > 6, trackParam.ellip<=3)
                warning('Requested ellipticity thresold is better to be close to 5 which means that at least 2 of the best focus plane should be consistent ([1 2 3 2 1])');
            end
            
            [traces,counter] = obj.trackParticles(trackParam);
            
            obj.traces.trace = traces;
            obj.traces.nTrace = counter;
            
        end
        
        function [traces] = get3DTraces(obj,calib,fitZParam,fittingType)
            
            %Extract 3D traces for each particles
            list = obj.particles.List;
            tracesIdx = obj.particles.traces;
            pxSize = obj.info.pxSize;
            zRange = fitZParam.zRange;
            %ellipt range used for fitting
            elliptRange = fitZParam.ellipRange(1):0.001:fitZParam.ellipRange(2);
            %we weigh the average later base on how much out of focus the
            %plane was.
            wRange1 = length(elliptRange(elliptRange<=1));
            wRange2 = length(elliptRange(elliptRange>=1));
            weight1 = linspace(1,5,wRange1);
            weight2 = linspace(5,1,wRange2);
            finalWeight = [weight1 weight2];
            
            %Calculation of x-y-z position occurs within the loop here
            traces = zeros(length(list),6,obj.particles.nTraces);
            for i = 1 : length(list)
                if ~isempty(list{i})
                    for j = 1 : length(list{i})
                        if ~isnan(tracesIdx{i}{j})
                            currentPart = list{i}{j};
                            %find indices to data in correct ellipticity range
                            id = and(currentPart.ellip>=elliptRange(1),...
                                currentPart.ellip<=elliptRange(end));
                            id2 = currentPart.fMetric > 75;
                            id = logical(id.*id2);
                            if all(id==0)
                                if i>3
                                    if and(traces(i-3,1,tracesIdx{i}{j})~=0,traces(i-2,1,tracesIdx{i}{j})~=0)
                                    else
                                        traces(i-3:i,:,tracesIdx{i}{j}) = 0;
                                    end
                                end
                                if i<length(list)-3
                                    if and(traces(i+3,1,tracesIdx{i}{j})~=0,traces(i+2,1,tracesIdx{i}{j})~=0)
                                    else
                                        traces(i:i+3,:,tracesIdx{i}{j}) = 0;
                                    end
                                end

                                
                            else
                                
                                ellip2Keep = currentPart.ellip(id);
                                idx = ellip2Keep;
                                for k = 1 :length(ellip2Keep)
                                    
                                    [~,idx(k)] = min(abs(elliptRange-ellip2Keep(k)));
                                    
                                    
                                end
                                
                                weight = finalWeight(idx);
                                %Weighed average
                                rowAvg = sum(diag(currentPart.row(id)* weight))/sum(weight) * pxSize;
                                colAvg = sum(diag(currentPart.col(id)* weight))/sum(weight) * pxSize;
                                %best focus Value
                                row = currentPart.row(3)* pxSize;
                                col = currentPart.col(3)* pxSize;
                                
                                if ~isempty(calib)
                                    %Calculation for Z is occur here
                                    [z,zAvg] = obj.getZPosition(currentPart,elliptRange,zRange,finalWeight,calib,fittingType);
                                    
                                else
                                    
                                    zAvg = 0;
                                    z = 0;
                                    
                                end
                                
                                if ~isnan(tracesIdx{i}{j})
                                    traces(i,:,tracesIdx{i}{j}) = [row col z rowAvg colAvg zAvg];
                                else
                                end
                                
                                
                            end
                        end
                    end
                end
            end
                      
        end
        
        function [zPos,zAvg] = getZPosition(obj,particle,elliptRange,zRange,finalWeight,calib,fittingType)
            %extract Z position based on ellipticity and zCalibrationration curve
            assert(~isempty(calib),'No zCalibration found, please give the zCalibration using giveZCal before calculating Z position');
            zCal = calib;
            relZ = obj.calibrated.oRelZPos;
            syncEllip = obj.zData.syncEllip;
            
            if ~isempty(syncEllip{particle.plane(3),1})
                zVec = zRange{particle.plane(3)}(1):zRange{particle.plane(3)}(2);
                switch fittingType
                    case 'poly'
                        fit = polyval(zCal{particle.plane(3),1},zVec);
                    case 'spline'
                        fit = ppval(zCal{particle.plane(3),2},zVec);%spline
                    otherwise
                        error('Unknown fitting type- only know "poly" and "spline"...')
                end

                %find the index of the value the closest to the particle
                %ellipticity
                [~,idx] = min(abs(fit-particle.ellip(3)));
                zPos = zVec(idx)+ relZ(particle.plane(3))*1000;

                id = and(particle.ellip>=elliptRange(1),...
                    particle.ellip<=elliptRange(end));
                data2Keep = particle(id,:);

                idx1 = zeros(size(data2Keep,1),1);
                zAvg = zeros(size(data2Keep,1),1);
                %Calculate z occur withing the loop
                for k = 1 :size(data2Keep,1)

                    [~,idx1(k)] = min(abs(elliptRange-data2Keep.ellip(k)));

                    zVectmp = min(zVec*2) : 1 : max(zVec*2);

                    switch fittingType
                    case 'poly'
                        fit = polyval(zCal{data2Keep.plane(k),1},zVectmp);
                    case 'spline'
                        fit = ppval(zCal{data2Keep.plane(k),2},zVectmp);%spline
                    otherwise
                        error('Unknown fitting type- only know "poly" and "spline"...')
                    end
                    
                    [~,idx] = min(abs(fit-data2Keep.ellip(k)));
                    zAvg(k) = zVectmp(idx) + relZ(data2Keep.plane(k))*1000;

                end

                weight = finalWeight(idx1);
                %weighed average occur here
                %diag is taken because of the matrix operation between weight
                %and zAvg
                [doAvg]  = obj.checkDoAverage;
                
                if doAvg
                    zAvg = sum(diag(zAvg(:)* weight(:)'))/sum(weight);
                else 
                    zAvg = zPos;
                end
            else
                zAvg = NaN;
                zPos = NaN;
            end
            
        end
        
    end
    
    methods (Static)
        
    end
    
    methods (Access = private)
        
    end
end