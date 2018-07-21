classdef zCalibration < handle
    %zCalibration is a class that will hold the information of several
    %zCalMovie and be able to process these information to create a
    %calibration curve together with displaying etc...
    
    properties
        path
        zCalMovies
        cal2D
        calib
        
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
        
        function retrieveZCalData(obj,detectParam)
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
                [traces,counter] = obj.zCalMovies.(['zCal' num2str(i)]).trackInZ;
                [zCalData] = obj.zCalMovies.(['zCal' num2str(i)]).getCalData(traces,counter);
                [~] = obj.zCalMovies.(['zCal' num2str(i)]).syncZCalData(zCalData);
                
            end
            disp('=================> DONE <===================');
            
        end
        
        function zCalibrate(obj)
            
            nfields = numel(fieldnames(obj.zCalMovies));
            allData = cell(8,3);
            %We regroup the data per plane of different movie together
            for i = 1: nfields
                
                for j = 1: length(obj.zCalMovies.(['zCal' num2str(i)]).zCalibration.syncEllip)
                    
                    allData{j,1} = [allData{j,1} ;obj.zCalMovies.(['zCal' num2str(i)]).zCalibration.syncEllip{j,1} ];
                    allData{1,2} = [allData{1,2} ;obj.zCalMovies.(['zCal' num2str(i)]).zCalibration.syncEllip{1,2} ];
                end
                
            end
            %Sort the data
            allData{1,3} = obj.zCalMovies.(['zCal' num2str(i)]).zCalibration.syncEllip{1,3};
            for j = 1: length(obj.zCalMovies.(['zCal' num2str(i)]).zCalibration.syncEllip)
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

