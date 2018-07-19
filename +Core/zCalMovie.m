classdef zCalMovie < Core.mpLocMovie
    %zCal will hold the information related to zzCalibrationration as well as all
    %the method linked to the zzCalibrationration.
    
    properties (SetAccess = 'private')
        zCalibration
        traces
    end
    
    methods
        
        function obj = zCalMovie(raw,cal,candidatePos)
            %UNTITLED6 Construct an instance of this class
            %   Detailed explanation goes here
            obj  = obj@Core.mpLocMovie(raw,cal);
            
            if nargin == 3
                
                    obj.candidatePos = candidatePos;
                    
            end
            
        end
        
        function set.zCalibration(obj, zCalibration)
            
            obj.zCalibration = zCalibration;
            
        end

        function [traces,counter] = trackInZ(obj)
            %track the particle in the Z direction (3rd dimension here)
            assert(~isempty(obj.calibrated),'Data should be calibrated to do ZzCalibrationration');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before zzCalibrationration');
            assert(~isempty(obj.particles), 'No particles found, please run superResConsolidate method before doing ZzCalibrationration');
            
            %We copy the List as boolean to keep track of where there are
            %still particles left
            [listBool] = obj.copyList(obj.particles.List,1);
            %We copy as NaN for storage of the traces;
            [traces]   = obj.copyList(obj.particles.List,NaN);
            %We pick the first particle available
            [idx] = obj.pickParticle(listBool);
            counter = 1;
            errCount =1;
            while (idx)
            %loop until there is no particle (pickParticle return false)   
            if errCount>1000
                warning('While loop ran for unexpectedly longer time');
                break;
                
            end
            %Connect particles (cf consolidation but across frames
            [listIdx] = obj.connectParticles(obj.particles.List,listBool,idx);
            %if the particle was connected in less than 5 frames we remove
            % its appearance from the list bool
            if length(listIdx) < 5 
                
                [listBool] = obj.removeParticles(listBool,listIdx);
                
            else
            %Otherwise we store traces, increment counter and remove.
            [traces]  = obj.storeTraces(traces,listIdx,counter);
            counter = counter +1;
            [listBool] = obj.removeParticles(listBool,listIdx);

            end
            % We pick a new particle and start all over again
            [idx] = obj.pickParticle(listBool);
            errCount = errCount +1;
            end
            counter = counter -1;
            obj.traces.trace = traces;
            obj.traces.nTrace = counter;
      
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
                            planes = obj.particles.List{i}{j}(:,end);
                            planes = planes (~isnan(planes));
                            
                            for k = 1 : length(planes)
                                data2Store = obj.particles.List{i}{j};
                                data2Store = data2Store(data2Store(:,end) == planes(k),:);
                                zCalData{planes(k),traces{i}{j}}(i,:) = [data2Store(3:end) i] ;
                                
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
                           
                            zCalData{i,j}(zCalData{i,j}(:,1)==0,:) = [];
                            
                       else
                       end
                       
                   end
               end
               obj.zCalibration.calData = zCalData;
        end 
        
        function [zData] = zCalibrate(obj)
            assert(~isempty(obj.calibrated),'Data should be calibrated to do ZzCalibrationration');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before zzCalibrationration');
            assert(~isempty(obj.particles), 'No particles found, please run superResConsolidate method before doing ZzCalibrationration');
            
            if ~isempty(obj.particles.Traces)
                quest = 'Some tracked traces were found in the object, do you want to keep them or run again ?';
                title = 'Question to User';
                btn1  = 'Keep';
                btn2 = 'Run again';
                defbtn = 'Keep';
                answer = questdlg(quest,title,btn1,btn2,defbtn);
                
                switch answer
                    case 'Keep'
                    case 'Run again'
                       
                        obj.trackInZ;
                        
                    otherwise
                        error('Unknown answer to question dialog ');
                        
                end
                
            else
                
                %We first track particle in Z
                obj.trackInZ;
                
            end
            %Transform the traces into zCalibrationration data (extracting data
            %per plane
            zCalData = obj.getCalData(obj.particles.Traces,obj.particles.nTraces);
            
            %synchronize the data (aligned them in Z)
            zSyncCalData = obj.syncZCalData(zCalData);
            %Calculate the zCalibrationration curve
            zCal =  Core.zCalMovie.calZCalibration(zSyncCalData);
            %store
            obj.zCalibration.cal = zCal;
            obj.zCalibration.calData = zCalData;
            obj.zCalibration.syncEllip = zSyncCalData;
            
            zData = obj.zCalibration;
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
        
        function showParticlesTracked(obj,ips)
            %display a movie where both localization and consolidation are
            %showed. The display needs to be improved !
            %TODO : improve display of the movie
            assert(~isempty(obj.particles.Traces),'You need to get the traces before displaying them');
            
            nPlanes = obj.calibrated.nPlanes;
            colors = rand(obj.particles.nTraces,3);
            
            for i = 1 : length(obj.particles.Traces)
            nParticles = length(obj.particles.Traces{i});
            obj.showCandidate(i);
            
            h = gcf;
            
                for j = 1 : nPlanes
                    subplot(2,nPlanes/2,j)
                    hold on
                    for k = 1 : nParticles
                        currPart = obj.particles.List{i}{k};
                        labelColor = obj.particles.Traces{i}{k};
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
        
        function showTraces(obj)
            %Display the x-y-z traces
            assert(~isempty(obj.particles.Traces),'You need to get the traces before displaying them');
            assert(isfield(obj.particles,'traces3D'),'You need to extract the 3D traces before display');
            traces = obj.particles.traces3D;
            npart = size(traces,3);
            %plot XYZ for every particles
            
            figure()
            
            for i = 1:npart
                data = traces(:,:,i);
                data = data(data(:,1)~=0,:);
                fprintf('std in X from best focus: %0.2f \n',nanmedian(nanstd(data(:,1))));
                fprintf('std in y from best focus: %0.2f \n',nanmedian(nanstd(data(:,2))));
                fprintf('std in Z from best focus: %0.2f \n',nanmedian(nanstd(data(:,3))));
                fprintf('std in X from mean planes: %0.2f \n',nanmedian(nanstd(data(:,4))));
                fprintf('std in y from mean planes: %0.2f \n',nanmedian(nanstd(data(:,5))));
                fprintf('std in Z from mean planes: %0.2f \n',nanmedian(nanstd(data(:,6))));
                
                
                subplot(2,1,1)
                hold on
                plot(data(:,1) - data(1,1));
                plot(data(:,2) - data(1,2));
                plot(data(:,3) - data(1,3));
                title('X Y Z position for different particle in best focus')
                xlabel('Frame')
                ylabel('Position(nm)')
                hold off
                
                subplot(2,1,2)
                hold on
                plot(data(:,4) - data(1,4));
                plot(data(:,5) - data(1,5));
                plot(data(:,6) - data(1,6));
                xlabel('Frame')
                ylabel('Position(nm)')
                title('X Y Z position for different particle mean')
                hold off
            end
               
            %plot Euclidian distance   
            figure()
            hold on
            for i = 1:npart
                data = traces(:,:,i);
                data = data(data(:,1)~=0,:);
                 %Calc euclidian distance
                eucl = sqrt((data(:,1)-data(1,1)).^2 + (data(:,2)-data(1,2)).^2 +...
                    (data(:,3)-data(1,3)).^2 );
                medEucl = sqrt((data(:,4)-data(1,4)).^2 + (data(:,5)-data(1,5)).^2 +...
                    (data(:,6)-data(1,6)).^2 );
                
                eucl2D = sqrt((data(:,1)-data(1,1)).^2 + (data(:,2)-data(1,2)).^2);
                medEucl2D = sqrt((data(:,4)-data(1,4)).^2 + (data(:,5)-data(1,5)).^2);
                
                fprintf('std in 2D from best focus: %0.2f \n',nanmedian(nanstd(eucl2D)));
                fprintf('std in 2D from mean of planes: %0.2f \n',nanmedian(nanstd(medEucl2D)));
                fprintf('std in 3D from best focus: %0.2f \n',nanmedian(nanstd(eucl)));
                fprintf('std in 3D from mean of planes: %0.2f \n',nanmedian(nanstd(medEucl)));
                
                subplot(2,2,1)
                hold on
                plot(eucl);
                title({'3D euclidian distance'; 'From best focus'})
                hold off
                
                subplot(2,2,2)
                hold on
                plot(medEucl);
                title({'3D euclidian distance'; 'From median'})
                hold off
                
                subplot(2,2,3)
                hold on
                plot(eucl2D);
                title({'2D euclidian distance'; 'From best focus'})
                hold off
                
                subplot(2,2,4)
                hold on
                plot(medEucl2D);
                title({'2D euclidian distance'; 'From median'})
                hold off
                
                %ylim([-400 400]);
                

            end
            hold off
            
            
        end
        
        function showZCalibration(obj)
            assert(~isempty(obj.zCalibration),'No zzCalibrationration found, please run z zCalibrationration before display');
            relZ = obj.calibrated.oRelZPos*1000;%in nm
            figure()
            for i = 1 : length(obj.zCalibration.syncEllip)
                
                z =  obj.zCalibration.syncEllip{i}(:,1)+relZ(i);
                ellip = obj.zCalibration.syncEllip{i}(:,2);
                
                ellip1 = ellip(ellip>=1);
                z1 = z(ellip>=1);
                
                ellip2 = 1./ellip(ellip<=1);
                z2 = z(ellip<=1);
               yAx = 1:0.1:2;
               subplot(2,1,1)
               hold on
               scatter(z1,ellip1)
               plot(ones(length(yAx),1)*relZ(i),yAx)
               title('Elliptiticy Elongated in Y')
               xlabel('zPos (nm)')
               ylabel('Ellipticity (sigY/sigX)')
               ylim([1 2])
               hold off
               
               subplot(2,1,2)
               hold on
               scatter(z2,ellip2)
               plot(ones(length(yAx))*relZ(i),yAx)
               title('Elliptiticy Elongated in X')
               xlabel('zPos (nm)')
               ylabel('Ellipticity (sigX/sigY)')
               ylim([1 2])
               hold off
                
            end
            
            figure()
            scatter(obj.zCalibration.syncEllip{1,2}(:,1),obj.zCalibration.syncEllip{1,2}(:,2));
            xlabel('ZPosition')
            ylabel('Ellipticity')
            title('Ellipticity-Z curve for all the planes superimposed')
            
            figure()
            
            for i = 1 : length(obj.zCalibration.syncEllip)
                
                dataCurrentPlane = obj.zCalibration.syncEllip{i};
                
                [binnedData] = Plotting.qBinning(dataCurrentPlane,length(dataCurrentPlane)/5);
                zVec = min(dataCurrentPlane(:,1)):max(dataCurrentPlane(:,1));
                %retrieving fit to display
                p = obj.zCalibration.cal{i};
                fit = polyval(p,zVec);
                %shifting according to the plane
                zVec = zVec + relZ(i) ;
                dataCurrentPlane(:,1) = dataCurrentPlane(:,1)+ relZ(i); 
                binnedData(:,1) = binnedData(:,1) +relZ(i);
                
                subplot(1,2,1)
                hold on
                scatter(binnedData(:,1),binnedData(:,2))
                plot(zVec,fit)
                title('Binned data fitted with polynomial')

                subplot(1,2,2)
                hold on
                scatter(dataCurrentPlane(:,1),dataCurrentPlane(:,2))
                plot(zVec,fit)
                title('Full data fitted with polynomial')
                
            end
            

        end
    
        function [zPos,zAvg] = getZPosition(obj,particle,elliptRange,finalWeight)
            %extract Z position based on ellipticity and zCalibrationration curve
            assert(~isempty(obj.zCalibration),'No zzCalibrationration found, please run z zCalibrationration calculating Z position');
            zCal = obj.zCalibration.cal;
            relZ = obj.calibrated.oRelZPos;
            syncEllip = obj.zCalibration.syncEllip;

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
        
        function [zSyncCalData] = syncZCalData(obj,zCalData)
            %Fit the ellipiticty - zPos data for each particles and
            %synchronized them by putting the point where ellipticity = 1
            %to zPos=0+PlanePosition. It uses the zCaldata which represent
            %the data point for each particle for a particle plane. Fitting
            %is done by polynomial to increase accuracy (if e = 1.05 then
            %this point will be locate not exactly at 0 but will be shifted
            deg = 4;
            minEllipt = 0.77;
            maxEllipt = 1.6;
            [zStep,~] = obj.getZPosMotor;
            zSyncCalData = cell(8,2);
%             figure
            for i = 1:size(zCalData,1)
                for j = 1:size(zCalData,2)
                    %Extract data in acceptable ellipticity range
                    %=>arbitrarily chosen when doing the fit, we stored it
                    %and reuse it here
                    if ~isempty(zCalData{i,j})
                        zPos = zCalData{i,j}(:,5)*zStep;
                        ellipt = zCalData{i,j}(:,1);
                        zPos = zPos(and(ellipt<maxEllipt, ellipt > minEllipt));
                        ellipt = ellipt(and(ellipt<maxEllipt, ellipt > minEllipt));

                        %Now we shift in z the value closest to ellipt =1 for
                        %the fit == rough synchronization    
                        [~,idx] = min(abs(ellipt-1));
                        focus1 = zPos(idx);
                        zPos = zPos - focus1;

                        if length(ellipt)<deg
                        else
                            %Do the fit and extract the exact z position of the focus
                            p = polyfit(zPos,ellipt,deg);
                            zVec = min(zPos):0.001:max(zPos);%1nm step
                            %this extraction has a 1nm accuracy
                            if or(min(zPos)>-0.1, max(zPos) <0.1)
                                zVec = -1:0.001:1; 
                            end

                            fit = polyval(p,zVec);
                            [~,idx] = min(abs(fit-1));
                            focus2 = zVec(idx);                  
                            shift = focus1+focus2;%take into account both synchronization
                            zCalData{i,j}(:,6) = ((zCalData{i,j}(:,5))*zStep -shift)*1000;%transform into nm keeping particle separated
                            fullSync = ((zCalData{i,j}(:,5))*zStep -shift)*1000;%mix all particle together
                            %tmp Store
                            zSyncCalData{i,1} = [zSyncCalData{i,1}; zCalData{i,j}(:,[6 1])];
                            zSyncCalData{1,2} = [zSyncCalData{1,2}; [fullSync zCalData{i,j}(:,1)]];
                        end
                    end
                end
                %tmp Store
                [~,ind] = sort(zSyncCalData{i,1}(:,1));
                zSyncCalData{i,1} = zSyncCalData{i,1}(ind,:);
                
            end
            %final storing and output
            [~,ind] = sort(zSyncCalData{1,2}(:,1));
            zSyncCalData{1,2} = zSyncCalData{1,2}(ind,:);
            zSyncCalData{1,3} = [minEllipt, maxEllipt, deg];
            obj.zCalibration.syncEllip = zSyncCalData;
        end
        
    end
    
    methods (Static)
        
        function [zCalibration] = calZCalibration(zSyncCalData)
            %function that take the synchronized z-ellip Data and fit, for
            %each planes with a polynomial. It stores the coeff of the
            %polynomials
            zCalibration = cell(length(zSyncCalData),1);
            deg = zSyncCalData{1,3}(3);
            
            for i = 1: length(zSyncCalData)
                dataCurrentPlane = zSyncCalData{i};
                [binnedData] = Plotting.qBinning(dataCurrentPlane,length(dataCurrentPlane)/5);
                
                zVec = min(dataCurrentPlane(:,1)):max(dataCurrentPlane(:,1));
     
                p = polyfit(dataCurrentPlane(:,1),dataCurrentPlane(:,2),deg);
                fit = polyval(p,zVec);
%                 
%                 figure
%                 subplot(1,2,1)
%                 hold on
%                 scatter(binnedData(:,1),binnedData(:,2))
%                 plot(zVec,fit)
%                 
%                 subplot(1,2,2)
%                 hold on
%                 scatter(dataCurrentPlane(:,1),dataCurrentPlane(:,2))
%                 plot(zVec,fit)
 
                zCalibration{i} = p;
                
            end
            
        end
        
        
        
    end
    

     methods (Access = private)
        
        
        
     end
end