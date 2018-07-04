classdef Movie < handle
    %General definition of a movie (master class from which many other
    %movie type will inherit
    
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
            
            assert(isfolder(raw), 'The given path is not a folder');
            %Check Given path
            [file2Analyze] = getOMETIF(obj,raw);
      
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
        
        function getRaw(obj,path)
            
            obj.raw = path;
            
        end
        
        function giveInfo(obj)
            
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
            
            assert(length(idx)==1,'Error too many frame requested, please load one at a time');
            
            [idx] = obj.checkFrame(idx);
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
            [idx] = obj.checkFrame(idx);
            %LoadCam
            [movC1,movC2,~] = Load.Movie.ome.load(obj.raw.frameInfo,obj.raw.movInfo,idx);
            data.Cam1 = movC1;
            data.Cam2 = movC2;
        end
        
        function [file2Analyze] = getFileInPath(~, path, ext)
           
            assert(ischar(path),'The given path should be a char');
            assert(ischar(ext),'The given extension should be a char');
            assert(isfolder(path),'The path given is not a folder')
            
            folderContent = dir(path);
            index2Images  = contains({folderContent.name},ext);
            file2Analyze  = folderContent(index2Images);
            
        end
        
        function [file2Analyze] = getOMETIF(obj,path)
            
            expExt = '.ome.tif';
            %Check Given path
            [file2Analyze] = obj.getFileInPath(path, expExt);
            assert(~isempty(file2Analyze),sprintf('The given directory does not any %s files',expExt));
            
        end
        
        function [frames]   = checkFrame(obj, frames) 
            
            testFrame = mod(frames,1);
            
            if all(testFrame<1e-4)
                
            else
                
                frames = round(frames);
                warning('Some frame were not integers, they were rounded');
                
            end
            
            assert(isvector(frames),'Frames should be a vector of integers');
            assert(max(frames) <= obj.raw.maxFrame(1),'Request exceeds max frame');
            assert(min(frames) >0, 'Indexing in matlab start from 1');
            
        end
        
    end
    
    methods (Access = private)
        
        
    end
    
end

