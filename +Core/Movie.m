classdef Movie < handle
    %General definition of a movie (master class from which many other
    %movie type will inherit
    %The Movie hold information and path to the movie but does not store
    %any data inside. Methods allow to display and do stuff with the data.
    
    properties (SetAccess = 'private')
       raw
        
    end
    properties 
       info 
    end
    
    methods
        %Constructor
        function obj = Movie(raw, info)
            %MOVIE Construct an instance of this class
            %   We allow the user to create this object with various number
            %   of input allowing therefore to restart the analysis at any
            %   steps in the process.
           
            %Give a status depending on input.
            switch nargin
                case 0 
                    
                    error('A path to the folder where the movie is is needed to create a movie object.')
                
                case 1 
                    
                     obj.raw = raw;
                     info.type = 'normal';
                     info.runMethod = 'load';
                     
                case 2
                    
                    if ~isfield(info,'runMethod')
                        quest = 'If we find data from previous analysis, do you want to load them or run the analysis again ?';
                        title = 'Question to User';
                        btn1  = 'Load';
                        btn2 = 'run again';
                        defbtn = 'Load';
                        answer = questdlg(quest,title,btn1,btn2,defbtn);
                        
                        switch answer
                            case 'Load'
                                
                                info.runMethod = 'load';
                                
                            case 'run again'
                                
                                info.runMethod = 'run';
                            otherwise
                                error('WTF');
                        end
                    end
                    
                    obj.raw = raw;
                    obj.info = info;
                                         
                otherwise 
                    
                    error('Unexpected number of argument');
                    
            end 
        end
        
        function set.raw(obj,raw)
            %This function will be adapted later to be able to take any
            %type of Movie (not only OME-TIFF).
            assert(isfolder(raw), 'The given path is not a folder');
            %Check Given path
            [file2Analyze] = Core.Movie.getOMETIF(raw);
      
            if length(file2Analyze)>1
                fprintf('More than one Tiff, loading only:\n %s', file2Analyze(1).name);
                fullPath = [file2Analyze(1).folder filesep file2Analyze(1).name];
                [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(fullPath);
                
                if iscell(frameInfo)
                    disp('Those tiff are multi-Images, we combine the info...')
                    [frameInfo, totFrame] = Load.Movie.ome.combineFrameInfo(frameInfo,false);
                    movInfo.indivFrame = movInfo.maxFrame;
                    movInfo.maxFrame = totFrame;
                end
            else
            
                fullPath = [file2Analyze(1).folder filesep file2Analyze(1).name];
                [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(fullPath);
                movInfo.indivFrame = movInfo.maxFrame;
                %Check info for 2 cam
               if length(movInfo.Cam) ~= 2
                   warning('Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
               end 
            end
            
            obj.raw.movInfo   = movInfo;
            obj.raw.frameInfo = frameInfo;
            obj.raw.fullPath  = fullPath;
            obj.raw.maxFrame  = movInfo.maxFrame;
        end
        
        function set.info(obj,inform)
            
            assert(isstruct(inform),'Information is expected to be a structure');
            names = fieldnames(inform);
            for i = 1:numel(fields(inform))

              obj.info.(names{i}) = inform.(names{i});

            end
        end
        
        function [raw] = getRaw(obj)
            
            raw = obj.raw;
            
        end
        
        function [info] = getInfo(obj)
            
            info = obj.info;
            
        end
        
        function giveInfo(obj)
            %Make a prompt asking some question to the user.
            prompt = {'Enter the pixel size: ','Enter the NA of the objective ',...
                'Enter the emission wavelength', 'Any comment about experiment?'};
            dlgTitle = 'Information about experimental parameters';
            numLines = 1;
            defaultVal = {'95','1.2','520',''};
            answer = inputdlg(prompt, dlgTitle,numLines,defaultVal);
            
            assert(~isempty(answer),'User canceled input dialog, Simulation was aborted')

            pxSize = str2double(answer(1));
            assert(~isnan(pxSize),'Number of Frame should be numerical');%If not a number
            
            NA = str2double(answer(2));
            assert(~isnan(NA),'NA should be numerical');
            
            emW = str2double(answer(3));
            assert(~isnan(emW),'Emission wavelength should be numerical');
            
            comment = answer(4);
            %Calculate some setup parameters
            sigma_nm = 0.25 * emW/NA;
            FWHMnm = sigma_nm * sqrt(8*log(2));
            FWHM_pix = FWHMnm/pxSize;
            sigmaPix = sigma_nm/pxSize;
            %store info        
            obj.info.pxSize = pxSize;
            obj.info.NA = NA;
            obj.info.emW = emW;
            obj.info.FWHM_px =  FWHM_pix;
            obj.info.sigma_px = sigmaPix;
            obj.info.comment = comment;

        end
        
        function h = showFrame(obj,idx,scaleBar)
            
            %To display a frame as a figure
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            pxSize = obj.info.pxSize/1000;%in µm
            scaleBarPx = scaleBar/pxSize;
            [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
            [frame] = getFrame(obj,idx);            
            assert(isstruct(frame),'Error unknown data format, data should be a struct');
            
            fNames = fieldnames(frame);
            idx2Empty = structfun(@isempty, frame);
            idx2Data = find(idx2Empty==0);
            nImages = length(find(idx2Empty));
                        
            h = figure(1);
            h.Name = sprintf('Frame %d',idx);
            
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
            
            for i = 1:nImages
                
                currentIM = frame.(fNames{idx2Data(i)});
                subplot(nImages,1,i)
                if strcmp(obj.info.type,'transmission')
                    colormap('gray');
                    %calculate reflectance (somehow better than absorbance)
                    %currentIM = imcomplement(currentIM);
                else
                    colormap('jet');
                end
                
                imagesc(currentIM)
                hold on
                x = size(currentIM,2)-scaleBarPx-(0.05*size(currentIM,2)):size(currentIM,2)-0.05*size(currentIM,2);
                y = ones(1,length(x))*size(currentIM,1)-0.05*size(currentIM,2);
                text(mean(x),mean(y)-0.05*size(currentIM,1),[num2str(scaleBar) ' µm'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',14);
                plot(x,y,'-w','LineWidth',5);

                axis image;
                caxis([min(min(min(currentIM))), max(max(max(currentIM)))]);
                
                a = gca;
                a.XTickLabel = [];
                a.YTickLabel = [];
                a.GridColor = [1 1 1];
                title({fNames{i}, sprintf('Frame %d',idx)});
                hold off
            end
        end
        
        function [data] = getFrame(obj,idx)
            
            assert(length(idx)==1,'Requested frame exceed the size of the movie');
             [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
            %LoadCam
            [movC1,movC2,~] = Load.Movie.ome.load(obj.raw.frameInfo,obj.raw.movInfo,idx);
            data.Cam1 = movC1;
            data.Cam2 = movC2;
        end
        
        function playMovie(obj)
            %TODO: Code a good way of playing the movie;
        end
        
        function saveMovie(obj,ext,frameRate,scaleBar,frames,plane)
            
            switch nargin
                case 4
                    nFrames = obj.raw.movInfo.maxFrame(1);
                    frames = 1:nFrames;
                    plane = [];
                case 5
                    
                    plane = [];
                case 6
                    
            end
            nFrames = length(frames);
            path2File = obj.raw.movInfo.Path;
            filename=sprintf('%s%sfullMovie.%s', path2File,'\',ext);

            for j = 1:nFrames
                if ~isempty(plane)
                    Fig = obj.showFrame(j,scaleBar,plane);
                else
                    Fig = obj.showFrame(j,scaleBar);
                end
                hold on
                %scale bar
                ax = gca;
                
                set(ax,'visible','off');
                set(gca,'position',[0 0 1 1],'units','normalized')
                axis image;
                drawnow;

                hold off
                frame = getframe(Fig);
                switch ext
                    case 'mp4'
                mov(j) = frame;
                    case 'gif'

                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);

                    if j == 1

                        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);

                    else

                        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');

                    end

                end
            end
            
            if strcmp(ext,'mp4')
                v = VideoWriter(filename,'MPEG-4');
                v.FrameRate = frameRate;
                open(v)
                writeVideo(v,mov);
                close(v)
            end
        end
        
    end
    
    methods (Static)
        
        function [frames]       = checkFrame(frames,maxFrame) 
            %Short method that make sure that the frame are making sense.
            testFrame = mod(frames,1);
            
            if all(testFrame<1e-4)
                
            else
                
                frames = round(frames);
                warning('Some frame were not integers, they were rounded');
                
            end
            
            assert(isvector(frames),'Frames should be a vector of integers');
            assert(max(frames) <= maxFrame(1),'Request exceeds max frame');
            assert(min(frames) >0, 'Indexing in matlab start from 1');
            
        end
        
        function [file2Analyze] = getFileInPath(path, ext)
            %Small method to extract the file of a certain extension in a
            %given path
            assert(ischar(path),'The given path should be a char');
            assert(ischar(ext),'The given extension should be a char');
            assert(isfolder(path),'The path given is not a folder')

            folderContent = dir(path);
            index2Images  = contains({folderContent.name},ext);
            file2Analyze  = folderContent(index2Images);
            
        end
        
        function [file2Analyze] = getOMETIF(path)
            %Variant of getFileInPath for .ome.tif file
            expExt = '.ome.tif';
            %Check Given path
            [file2Analyze] = Core.Movie.getFileInPath(path, expExt);
            assert(~isempty(file2Analyze),sprintf('The given directory does not any %s files',expExt));
            
        end
        
    end
    
    methods (Access = private)
        
        
        
    end
    
end

