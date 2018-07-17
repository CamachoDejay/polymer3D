classdef Movie < handle
    %General definition of a movie (master class from which many other
    %movie type will inherit
    %The Movie hold information and path to the movie but does not store
    %any data inside. Methods allow to display and do stuff with the data.
    
    properties (SetAccess = 'private')
        raw
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
                     
                case 2
                    
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
                
                fprintf('More than one Tiff, Loading %s \n', file2Analyze(1).name);
                
            end
            
            fullPath = [file2Analyze.folder filesep file2Analyze(1).name];
            [frameInfo, movInfo, ~ ] = Load.Movie.ome.getInfo(fullPath);
            %Check info for 2 cam
            assert(length(movInfo.Cam) == 2,'Only 1 camera found in the selected file, code only works with 2 cameras, will be updated later.');
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
            inform.pxSize = pxSize;
            inform.NA = NA;
            inform.emW = emW;
            inform.FWHM_px =  FWHM_pix;
            inform.sigma_px = sigmaPix;
            inform.comment = comment;
            
            obj.info = inform;
            
        end
        
        function showFrame(obj,idx)
            %To display a frame as a figure
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            
            [idx] = Core.Movie.checkFrame(idx,obj.raw.maxFrame(1));
            [frame] = getFrame(obj,idx);            
            assert(isstruct(frame),'Error unknown data format, data should be a struct');
            
            nImages = numel(fields(frame));
            fNames = fieldnames(frame);
            nsFig = 2;                
            h = figure(1);
            h.Name = sprintf('Frame %d',idx);
            
            for i = 1:nImages
                
                subplot(2,nImages/nsFig,i)
                imagesc(frame.(fNames{i}))
                axis image;
                grid on;
                a = gca;
                a.XTickLabel = [];
                a.YTickLabel = [];
                a.GridColor = [1 1 1];
                title({fNames{i}, sprintf('Frame %d',idx)});
                colormap('jet')
                
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

