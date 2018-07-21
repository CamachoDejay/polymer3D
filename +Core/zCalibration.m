classdef zCalibration < handle
    %zCalibration is a class that will hold the information of several
    %zCalMovie and be able to process these information to create a
    %calibration curve together with displaying etc...
    
    properties
        path
        zCalMovies
        cal2D
        calib
        traces3D
        zPosMotor
    end
    
    methods
        
        function obj = zCalibration(path2zCal,cal2D)
            %zCalibration Construct an instance of this class
            %   Detailed explanation goes here
            obj.path = path2zCal;
            obj.cal2D = cal2D;
        end
        
        function retrieveZCalMov(obj)
            %we get the zCalibration directory
            folder2Mov = dir(obj.path);
            %loop through the content of the directory
            for i = 3:size(folder2Mov,1)
                %If element i is a folder
                if folder2Mov(i).isdir
                    %Check if the directory
                    currDir = dir([folder2Mov(i).folder filesep folder2Mov(i).name]);
                    idx = contains({currDir.name}, 'ome.tif');
                    if ~all(idx==0)
                        obj.zCalMovies.(['zCal' num2str(i-2)]) = Core.zCalMovie([folder2Mov(i).folder filesep folder2Mov(i).name], obj.cal2D);
                    else
                        warning([folder2Mov(i).folder filesep folder2Mov(i).name ' did not contain any ome.Tif and is therefore ignored']);
                    end
                end
            end
        end
        
        function retrieveZCalData(obj,detectParam, fitZParam, trackParam)
            nfields = numel(fieldnames(obj.zCalMovies));
            for i = 1: nfields
                disp(['Retrieving data from zCal file ' num2str(i) ' / ' num2str(nfields) ' ...']);
                if i == 1
                    
                    obj.zCalMovies.(['zCal' num2str(i)]).giveInfo;
                    
                else
                    
                    obj.zCalMovies.(['zCal' num2str(i)]).info = obj.zCalMovies.(['zCal' num2str(1)]).getInfo;
                    
                end
                
                obj.zCalMovies.(['zCal' num2str(i)]).findCandidatePos(detectParam);
                obj.zCalMovies.(['zCal' num2str(i)]).superResConsolidate;
                [traces,counter] = obj.zCalMovies.(['zCal' num2str(i)]).trackInZ(trackParam);
                [zCalData] = obj.zCalMovies.(['zCal' num2str(i)]).getCalData(traces,counter);
                [~] = obj.zCalMovies.(['zCal' num2str(i)]).syncZCalData(zCalData,fitZParam);
                
            end
            disp('=================> DONE <===================');
            obj.calib.fitZParam = fitZParam;
            obj.calib.trackParam = trackParam;
        end
        
        function zCalibrate(obj)
            
            nfields = numel(fieldnames(obj.zCalMovies));
            allData = cell(8,3);
            %We regroup the data per plane of different movie together
            for i = 1: nfields
                
                for j = 1: length(obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip)
                    
                    allData{j,1} = [allData{j,1} ;obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip{j,1} ];
                    
                end
                
                allData{1,2} = [allData{1,2} ;obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip{1,2} ];
            end
            %Sort the data
            allData{1,3} = obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip{1,3};
            for j = 1: length(obj.zCalMovies.(['zCal' num2str(i)]).zData.syncEllip)
                [~,ind] = sort(allData{j,1}(:,1));
                allData{j,1} = allData{j,1}(ind,:);
                
            end
            %Sort all data together
            [~,ind] = sort(allData{1,2}(:,1));
            allData{1,2} = allData{1,2}(ind,:);
            %Fitting of the data occurs here
            zCalibration = obj.calZCalibration(allData);
            %storage to the object
            obj.calib.file = zCalibration;
            obj.calib.data = allData;
            
        end
        
        function showZCalibration(obj)
            
            assert(~isempty(obj.calib),'No zzCalibrationration found, please run z zCalibrationration before display');
            relZ = obj.zCalMovies.('zCal1').calibrated.oRelZPos*1000;%in nm
            figure()
            for i = 1 : length(obj.calib.data)
                
                z =  obj.calib.data{i}(:,1)+relZ(i);
                ellip = obj.calib.data{i}(:,2);
                
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
            scatter(obj.calib.data{1,2}(:,1),obj.calib.data{1,2}(:,2));
            xlabel('ZPosition')
            ylabel('Ellipticity')
            title('Ellipticity-Z curve for all the planes superimposed')
            
            figure()
            
            for i = 1 : length(obj.calib.data)
                
                dataCurrentPlane = obj.calib.data{i};
                
                [binnedData] = Plotting.qBinning(dataCurrentPlane,length(dataCurrentPlane)/5);
                zVec = min(dataCurrentPlane(:,1)):max(dataCurrentPlane(:,1));
                %retrieving fit to display
                p = obj.calib.file{i};
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
        
        function evalAccuracy(obj)
            
            nfields = numel(fieldnames(obj.zCalMovies));
            zMotor = cell(nfields,1);
            trace3D  = cell(nfields,1);
            
            for i = 1: nfields
                
                disp(['Retrieving 3D traces ' num2str(i) ' / ' num2str(nfields) ' ...']);
                [~,zMotor{i}] = obj.zCalMovies.(['zCal' num2str(i)]).getZPosMotor;
                trace3D{i} = obj.zCalMovies.(['zCal' num2str(i)]).get3DTraces (obj.calib.file,obj.calib.fitZParam);
                
            end
            obj.traces3D = trace3D;
            obj.zPosMotor = zMotor;
            disp('=================> DONE <===================');
        end
        
        function showAccuracy(obj)
       
            %Display the x-y-z traces
            assert(~isempty(obj.traces3D),'You need to get the traces before displaying them, use evalAccuracy to get the traces');
            trace = obj.traces3D;
            motor = obj.zPosMotor;
            accuracyFocus = [];
            accuracyMean = [];
            %plot XYZ for every particles       
            nfields = numel(fieldnames(obj.zCalMovies));
            figure()
            for i = 1: nfields
                
                currentTrace = trace{i};
                currentMotor = motor{i};
                currentMotor = (currentMotor - currentMotor(1))*1000;
                npart = size(currentTrace,3);
                
                for j = 1:npart
                    data = currentTrace(:,:,j);
                    if ~all(data==0)

                        data = data(data(:,1)~=0,:);
                        accuracyF2plot = (data(:,3) - data(find(data(:,3),1,'first'),3)) - currentMotor(1:size(data,1));
                        accuracyM2plot  = (data(:,6) - data(find(data(:,3),1,'first'),6)) - currentMotor(1:size(data,1));
                        
                        accuracyMean = [accuracyMean std(accuracyM2plot)];
                        accuracyFocus = [accuracyFocus std(accuracyF2plot)];
                          
                        subplot(2,2,1)
                        hold on
                        scatter(1:size(data(:,3),1),data(:,3) - data(1,3));
                        plot(currentMotor,'-b');

                        title('Z position for different particle in best focus')
                        xlabel('Frame')
                        ylabel('Position(nm)')
                        hold off

                        subplot(2,2,2)
                        hold on
                        scatter(1:size(data,1),data(:,6) - data(1,6));
                        plot(currentMotor,'-b');
                        xlabel('Frame')
                        ylabel('Position(nm)')
                        title('Z position for different particle mean')
                        hold off

                        subplot(2,2,3)
                        hold on
                        plot(nonzeros(accuracyF2plot))
                        xlabel('Frame')
                        ylabel('tracked Distance - motor movement')
                        title(' Z position for different particle in best focus')
                        hold off

                        subplot(2,2,4)
                        hold on
                        plot(nonzeros(accuracyM2plot))
                        xlabel('Frame')
                        ylabel('tracked Distance - motor movement')
                        title('Z position for different particle mean')
                        hold off
                    end
                end
            end
            
            disp(['Accuracy in Z for best focus is ' num2str(mean(accuracyFocus))]);
            disp(['Accuracy in Z for mean  ' num2str(mean(accuracyMean))]);
            
%             %plot Euclidian distance   
%             figure()
%             hold on
%             for i = 1:npart
%                 data = trace(:,:,i);
%                 data = data(data(:,1)~=0,:);
%                  %Calc euclidian distance
%                 eucl = sqrt((data(:,1)-data(1,1)).^2 + (data(:,2)-data(1,2)).^2 +...
%                     (data(:,3)-data(1,3)).^2 );
%                 medEucl = sqrt((data(:,4)-data(1,4)).^2 + (data(:,5)-data(1,5)).^2 +...
%                     (data(:,6)-data(1,6)).^2 );
%                 
%                 eucl2D = sqrt((data(:,1)-data(1,1)).^2 + (data(:,2)-data(1,2)).^2);
%                 medEucl2D = sqrt((data(:,4)-data(1,4)).^2 + (data(:,5)-data(1,5)).^2);
%                 
%                 fprintf('std in 2D from best focus: %0.2f \n',nanmedian(nanstd(eucl2D)));
%                 fprintf('std in 2D from mean of planes: %0.2f \n',nanmedian(nanstd(medEucl2D)));
%                 fprintf('std in 3D from best focus: %0.2f \n',nanmedian(nanstd(eucl)));
%                 fprintf('std in 3D from mean of planes: %0.2f \n',nanmedian(nanstd(medEucl)));
%                 
%                 subplot(2,2,1)
%                 hold on
%                 plot(eucl);
%                 title({'3D euclidian distance'; 'From best focus'})
%                 hold off
%                 
%                 subplot(2,2,2)
%                 hold on
%                 plot(medEucl);
%                 title({'3D euclidian distance'; 'From median'})
%                 hold off
%                 
%                 subplot(2,2,3)
%                 hold on
%                 plot(eucl2D);
%                 title({'2D euclidian distance'; 'From best focus'})
%                 hold off
%                 
%                 subplot(2,2,4)
%                 hold on
%                 plot(medEucl2D);
%                 title({'2D euclidian distance'; 'From median'})
%                 hold off
%                 
%                 %ylim([-400 400]);
%                 
% 
%             end
%             hold off
                     
        end
    end
    
    methods (Access = private)
        
        function [zCalibration] = calZCalibration(~,zSyncCalData)
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
                
                zCalibration{i} = p;
                
            end
            
        end
        
    end
end

