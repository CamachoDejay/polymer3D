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
        
        function [traces] = getTraces(obj)
            assert(~isempty(obj.traces3D),'please run the tracking before getting the traces');
            traces = obj.traces3D;
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
                if errCount>10000
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
            
            filename =[obj.raw.movInfo.Path filesep 'Traces3D.mat'];
            
            save(filename,'trace3D');
            
            
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
        
        function evalAccuracy(obj,dim)
            [xStep,xMotor] = obj.getXPosMotor;
            [yStep,yMotor] = obj.getYPosMotor;
            [zStep,zMotor] = obj.getZPosMotor;
            
            switch nargin
                case 1
                    
                    dim = {'row','col','z'};
                    obj.getAccuracy(xMotor,dim{2});
                    obj.getAccuracy(yMotor,dim{1});
                    obj.getAccuracy(zMotor,dim{3});
         
                case 2
                    
                    switch dim
                        case 'x'
                            dim = 'col';
                            obj.getAccuracy(xMotor,dim{2});
                        case 'y'
                            dim = 'row';
                            obj.getAccuracy(yMotor,dim{1});
                        case 'z'
                            obj.getAccuracy(zMotor,dim{3});
                    end
                    
                otherwise
                    error('too many input arguments');
            end
            
            if all([xStep, yStep, zStep] ==0)
                warning(['no movement of the motor to compare with,',...
                    'we will assume the sample does not move for evaluating accuracies',...
                    'if it is not the case, it might result in unexpected numbers']);
            end           
        end
        
        function [int,SNR] = getAvgIntensity(obj)
            assert(~isempty(obj.traces3D),'You need to extract 3D traces before extracting average intensity');
            traces = obj.traces3D;
            nTraces = length(traces);
            int = zeros(nTraces,1);
            SNR = zeros(nTraces,1);
            for i = 1: length(traces)
                currentTrace = traces{i};
                int(i) = mean(currentTrace.intensity);
                SNR(i) = mean(currentTrace.SNR);
                
            end
            
            int = mean(int);
            SNR = mean(SNR);
           
            
        end
        
        function getTracesMovie(obj,frames,idx2Trace,ROI,frameRate)
            assert(~isempty(obj.traces3D),'You need to extract 3D traces before extracting movies');
            obj.getPartMovie(obj,frames,idx2Trace,ROI,frameRate);
            
            obj.getTrace3DMovie(obj,frames,idx2Trace,frameRate);
        end
        
        function getPartMovie(obj,frames,idx2Trace,ROI,frameRate,scaleBar)
        assert(~isempty(obj.traces3D),'You need to extract 3D traces before extracting particle trace movie');
        if nargin <6
            scaleBar = 500;
        elseif nargin <5
            error('not enough input arguments');
        else
            error('too many input arguments');
        end
        path2File = obj.raw.movInfo.Path;
        traces = obj.traces3D;
        roiRadius = ROI;
        pxSize = obj.info.pxSize;
        currentTraces = traces {idx2Trace};
        mainPos = [round(mean(currentTraces.row)/pxSize) round(mean(currentTraces.col(1)/pxSize))];
        nFrames = length(frames);
        frames = currentTraces.frame(1:nFrames);
        scaleBarPx = scaleBar/pxSize;
        pos.row = currentTraces.row/pxSize - mainPos(1) + roiRadius + 1;
        pos.col = currentTraces.col/pxSize - mainPos(2) + roiRadius + 1;
        
        for i = 1:obj.calibrated.nPlanes

            currentPlane = obj.getPlane(i);
            % ROI = currentPlane;
            ROI = currentPlane(mainPos(1)-roiRadius:mainPos(1)+roiRadius,...
            mainPos(2)-roiRadius:mainPos(2)+roiRadius,:);
            mov = struct('cdata', cell(1,nFrames), 'colormap', cell(1,nFrames));
            Fig = figure;
            %to get as less white border as possible
            ax = gca;
            
            outerpos = ax.OuterPosition;
            ti = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
            
            for j = 1:nFrames
            
            gcf;
           
            imagesc(ROI(:,:,frames(j)))
            hold on
            %scale bar
            x = size(ROI,2)-scaleBarPx-(0.05*size(ROI,2)):size(ROI,2)-0.05*size(ROI,2);
            y = ones(1,length(x))*size(ROI,1)-0.05*size(ROI,2);
            text(mean(x)-2,mean(y)-2,[num2str(scaleBar) ' nm'],'Color','white','fontWeight','bold');
            plot(x,y,'-w','LineWidth',5);
            plot(pos.col(1:j),pos.row(1:j),'-b')
            caxis([min(min(min(ROI))) max(max(max(ROI)))]);
            set(ax,'visible','off');
            axis image;
            
            colormap('hot')
            drawnow;
            
            hold off
            mov(j) = getframe(Fig);

            end
            ext='.mp4';
            filename=sprintf('%s%splane%d-Trace%d%s', path2File,'\',i,idx2Trace,ext);
            v = VideoWriter(filename,'MPEG-4');
            v.FrameRate = frameRate;
            open(v)
            writeVideo(v,mov);
            close(v)

        end

        end

        function getTraces3DMovie(obj,frames,idx2Trace,frameRate)
        assert(~isempty(obj.traces3D),'You need to extract 3D traces before getting traces Movie');
        %             sizeMarker = 5;
        Fig = figure;

        path2File = obj.raw.movInfo.Path;
        traces = obj.traces3D;
        currentTraces = traces {idx2Trace};

        nFrames = length(frames);
        [frames] = obj.checkFrame(frames,size(currentTraces,1));
        frames = currentTraces.frame(1:nFrames);

        xAx = [min(currentTraces.col-mean(currentTraces.col)),...
            max(currentTraces.col-mean(currentTraces.col))];
        yAx = [min(currentTraces.row-mean(currentTraces.row)),...
            max(currentTraces.row-mean(currentTraces.row))];
        zAx = [min(currentTraces.z-mean(currentTraces.z)),...
            max(currentTraces.z-mean(currentTraces.z))];
        mov = struct('cdata', cell(1,nFrames), 'colormap', cell(1,nFrames));
        gcf;
        hold on
        for j = 1:nFrames

            colPlot = currentTraces.col(1:j) - mean(currentTraces.col);
            rowPlot = currentTraces.row(1:j) - mean(currentTraces.row);
            zPlot = currentTraces.z(1:j) - mean(currentTraces.z);
            %plotting with z coloring:
            patch([colPlot nan(size(colPlot))],[rowPlot nan(size(colPlot))],...
                [zPlot nan(size(colPlot))],[zPlot nan(size(colPlot))],...
                'EdgeColor','interp','FaceColor','none')

            xlim(xAx)
            ylim(yAx)
            zlim(zAx)
            view(3);

            mov(j) = getframe(Fig);

            xlabel('x Position (nm)');
            ylabel('y Position(nm)');
            zlabel('z Position(nm)');
        end

        ext='.mp4';
        filename=sprintf('%s%sTracking-Trace%d%s', path2File,'\',idx2Trace,ext);
        v = VideoWriter(filename,'MPEG-4');
        v.FrameRate = frameRate;
        open(v)
        writeVideo(v,mov);
        close(v)
        end
        
        function [RMSD,traces] = getRMSD(obj,dimension)
            
            assert(~isempty(obj.traces3D),'You need to extract 3D traces before extracting RMSD');
            
            traces = obj.traces3D;
            
            switch nargin
                case 1
                    
                    dimension = '3D';
                
                case 2
                    
                otherwise
                    
                    error('too many input arguments');
            
            end
            RMSD = cell(traces);
            
            for i = 1 : length(traces)
                
                currentTrace = traces{i};
                coord = [currentTrace.col,currentTrace.row,currentTrace.z];
                [MSD,~] = Core.MPTrackingMovie.calcMeanSqrD(coord,dimension);
                currentTrace.RMSD = zeros(size(currentTrace.row,1),1);
                currentTrace.RMSD(1:end-1) = MSD;
                traces{i} = currentTrace;
                RMSD{i} = MSD;
            end
            
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
        
        function [RMSD,D] = calcMeanSqrD(coord,dimension)
            %[D MSD]=meansqrD(cod, stp)
        % the function gives the distance ('D') between each point and
        % mean-squre-dispracement ('MSD')
        %
        %input parameters; cod is 2D or 3D aray including coordinate of trajectry.
        if size(coord,1)< size(coord,2)

            coord = coord';
            warning('vector was wrongly oriented so we took the transpose');

        end

        dim = size(coord,2);

        switch dim
            case 1

                coord(:,2:3) = 0;

            case 2

                coord(:,3) = 0;

            case 3


            otherwise

                error('unexpected dimension for the vector')

        end
        
        if strcmp(dimension,'2D')
            
            coord(:,3) = 0;
            
        elseif strcmp(dimension,'3D')
            
            assert(size(coord,2) ==3,'error with the dimension of the vector, expected 3D vector');
            
        end

        %%%Fist calculate distance between each frame

        DX = diff(coord(:,1).');
        DY = diff(coord(:,2).');
        DZ = diff(coord(:,3).');
        D = sqrt(DX.^2 + DY.^2 + DZ.^2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        RMSD = zeros(length(coord)-1,1);
        %Calculate mean-squere-displacement
        for i = 1:length(coord)-1

            stp = i;
            cnt =  1;
            D1  = [];
            while cnt<=stp && cnt+stp<=length(coord)

                idx = cnt:stp:length(coord);
                DX  = diff(coord(idx,1).');
                DY  = diff(coord(idx,2).');
                DZ  = diff(coord(idx,3).');
                D1  = [D1 sqrt(DX.^2 + DY.^2 + DZ.^2)];
                cnt = cnt+1;

                if ~isempty(D1)

                    D2=D1(~isnan(D1));

                    if ~isempty(D2)

                        RMSD(i) = mean(D2.^2);

                    else

                        RMSD(i) = NaN;

                    end
                end
            end %while
        end

        D   = D(:);
        RMSD = RMSD(:);
            
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
                            {'row','col','z','intensity','SNR','frame'}) =...
                            [currentPart(:,{'row','col','z','intensity','SNR'}),  table(i)];
                        
                        
                        fCounter(currentTrace) = fCounter(currentTrace)+1;
                    else
                    end
                    
                
                end
            end
            
            
            
        end
        
        function getAccuracy(obj,Motor,dim)
            assert(ischar(dim),'dimension is expected to be a char (row, col or z')
            traces = obj.traces3D;
            
            meanErr = zeros(length(traces),1);
            absMeanErr = meanErr;
            stdErr  = meanErr;
            Fig = figure;
            allTraces = zeros(length(Motor),length(traces));
            for i = 1:length(traces)
                
                currentTrace = traces{i};
                motorPos = Motor(currentTrace.frame(1:end))*1000;
                                                
                mot = motorPos - mean(motorPos);
                %reflect y axis because of inverted directions
                if or(strcmp(dim,'col'),strcmp(dim,'row'))
                   mot = motorPos - 2*(motorPos.*(mean(motorPos)/max(motorPos)))*mean(motorPos)/max(motorPos);
                   mot = mot - mean(mot);
                end
                
                traceErr = currentTrace.(dim) - mean(currentTrace.(dim));
                meanErr(i)    = mean(traceErr - mot);
                absMeanErr(i) = mean(abs(traceErr - mot));
                stdErr(i)     = std(traceErr - mot);
                allTraces(currentTrace.frame(1:end),i) = traceErr;
                if size(currentTrace,1) > length(Motor)/2
                    subplot(1,2,1)
                    hold on
                    scatter(currentTrace.frame,traceErr)
                    plot(currentTrace.frame,mot,'-r','Linewidth',2)
                    xlabel('frame')
                    ylabel([dim,' pos (nm) '])
                    title([dim ' localization compared with motor']);
                    
                    subplot(1,2,2)
                    hold on
                    plot(traceErr-mot);
                    
                    xlabel('frame')
                    ylabel([dim,' error compare to motor (nm) '])
                    title(['tracking error']);
                    hold off
                end

            end
            disp(['mean accuracy ', dim,': ', num2str(mean(meanErr)), ' nm']);
            disp(['abs mean ', dim,': ', num2str(mean(absMeanErr)), ' nm']);
            disp(['std ',dim,': ', num2str(mean(stdErr)), ' nm']);
            filename = [obj.raw.movInfo.Path filesep dim '-Fig'];
            saveas(Fig,filename,'svg');
            
            test = allTraces~=0;
            
            test = sum(test,1);
            [~,idx] = max(test);
            idx2Use = test==test(idx);
            allTraces = allTraces(:,idx2Use);
            Motor = Motor(allTraces(:,1)~=0);
            allTraces = allTraces(allTraces(:,1)~=0,:);
            
            Fig = figure;
            x = 1:length(Motor);
            
            
            if or(strcmp(dim,'col'),strcmp(dim,'row'))
                   mot = Motor - 2*(Motor.*(mean(Motor)/max(Motor)))*mean(Motor)/max(Motor);
            else
                mot = Motor;
            end
            
            y = (mot-mean(mot))*1000 ;
            err = mean(abs(allTraces-y),2);
            
            H   = Plotting.shadedErrorBar(x(:),y(:),err(:));
            
            xlabel('Frame');
            ylabel([dim,'Tracking position']);
            filename = [obj.raw.movInfo.Path filesep dim '-ShadedErrFig'];
            saveas(Fig,filename,'svg');
        end
        

    end
end

