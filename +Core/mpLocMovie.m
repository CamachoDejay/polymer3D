classdef mpLocMovie < Core.mpMovie
    %mpLocMovie hold all the method and info necessary for localization,
    %expension of this object will be first centered around Tracking but
    %it will also be easy to extend to STORM/PALM.
    
    properties
        
        candidatePos
        
    end
    
    methods
        function obj = mpLocMovie(raw, cal)
            
            obj  = obj@Core.mpMovie(raw,cal);
         
        end
        
        function set.candidatePos(obj,candidatePos)
            
            assert(iscell(candidatePos), 'Expecting a cell containing candidate positions');
            obj.candidatePos = candidatePos;
            
        end
        
        function findCandidatePos(obj,detectParam, frames)
            
            assert(~isempty(obj.info), 'Missing information about setup to be able to find candidates, please use giveInfo method first');
            assert(nargin>1,'not enough input argument');
            [file2Analyze] = getFileInPath(obj, obj.raw.movInfo.Path, '.mat');
            %Check if some candidate exists already
            [run, candidate] = obj.existCandidate(file2Analyze);
            
            if run
                switch nargin 

                    case 2

                        frames = 1: obj.calibrated.nFrames;
                        run = true;
                        disp('Running detection on every frame');

                    case 3

                        run= true;
                        [frames] = obj.checkFrame(frames);

                    otherwise

                        error('too many inputs');

                end

                [candidate] = obj.detectCandidate(detectParam,frames);
                       
            elseif ~isempty(candidate)
            else

                disp('getCandidatePos is a function that detects features in a movie');
                disp('To work, it needs to receive a structure containing 2 detection parameter:');
                disp('delta which is the radius around which detection is performed usually 6 pixels');
                disp('chi2 which characterize how certain you want to be about the fact that there is a molecule there');
                disp('Typically between 20 and 200');

            end
          
            obj.candidatePos = candidate;
            
        end
        
        function [candidate] = getCandidatePos(obj, frames)
            
            [idx] = obj.checkFrame(frames);
            candidate = obj.candidatePos{idx};
            
            if isempty(candidate)
                
                warning('There was no candidate found in this frames, please check that you ran findCandidate on this frame, if yes, check the data');
            
            end
        end
        
        function showCandidate(obj,idx)
            
            assert(length(idx)==1, 'Only one frame can be displayed at once');
            
            [idx] = obj.checkFrame(idx);
            assert(~isempty(obj.candidatePos{idx}),'There is no candidate found in that frame, check that you ran the detection for that frame');
            
            [frame] = getFrame(obj,idx);
            assert(isstruct(frame),'Error unknown data format, data should be a struct');
            
            nImages = numel(fields(frame));
            fNames = fieldnames(frame);
            nsFig = 2;
            
            candidate = obj.getCandidatePos(idx);
            rowPos    = candidate(:,1);
            colPos    = candidate(:,2);
            planeIdx  = candidate(:,3);
            
            h = figure(2);
            h.Name = sprintf('Frame %d',idx);
            for i = 1:nImages
                
                subplot(2,nImages/nsFig,i)
                hold on
                imagesc(frame.(fNames{i}))
                plot(colPos(planeIdx==i),rowPos(planeIdx==i),'g+')
                axis image;
                grid on;
                a = gca;
                a.XTickLabel = [];
                a.YTickLabel = [];
                a.GridColor = [1 1 1];
                title({fNames{i},sprintf(' Zpos = %0.3f',obj.calibrated.oRelZPos(i))});
                colormap('jet')
                hold off
                
            end
        end
        
    end
    
    methods (Access = private)
        
        function [run,candidate] = existCandidate(obj,file2Analyze)
            
            if any(contains({file2Analyze.name},'candidatePos')==true)
            quest = 'Some candidate were found in the raw folder, do you want to load them or run again ?';
            title = 'Question to User';
            btn1  = 'Load';
            btn2 = 'run again';
            defbtn = 'Load';
            answer = questdlg(quest,title,btn1,btn2,defbtn);
            
            switch answer
            case 'Load'
                
                candidate = load([file2Analyze(1).folder filesep 'candidatePos.mat']);
                candidate = candidate.candidate;
                run = false;
                
            case 'run again'
                
                run = true;
                candidate =[];
                
            end
            
            else
                
                run = true;
                candidate =[];
            end    
        end
        
        function [candidate] = detectCandidate(obj,detectParam,frames)
            
            assert(~isempty(obj.calibrated),'Data should be calibrated to detect candidate');
            assert(isstruct(detectParam),'Detection parameter should be a struct with two fields');
            nFrames = length(frames);
            currentCandidate = obj.candidatePos;

            if(isempty(currentCandidate))

                candidate = cell(obj.calibrated.nFrames,1);

            else

                candidate = currentCandidate;

            end

            %parameter for localization
            FWHM_pix = obj.info.FWHM_px;
            delta  = detectParam.delta;
            chi2   = detectParam.chi2;
            h = waitbar(0,'detection of candidates...');

            for i = 1 : 1:nFrames

                position = zeros(200,3);
                [volIm] = obj.getFrame(frames(i)); 
                nameFields = fieldnames(volIm);

                for j = 1:length(nameFields)

                    [ pos, ~, ~ ] = Localization.smDetection( double(volIm.(nameFields{j})),...
                        delta, FWHM_pix, chi2 );
                    startIdx = find(position==0,1,'First');
                    pos(:,3) = j;
                    position(startIdx:startIdx+size(pos,1)-1,:) = pos;

                end

                idx = find(position==0,1,'First');
                candidate{frames(i)} = position(1:idx-1,:);
                waitbar(i/nFrames,h,...
                    sprintf('detection of candidates in Frame %d/%d done',i,nFrames));
            end

            close(h);

            fileName = sprintf('%s%scandidatePos.mat',obj.raw.movInfo.Path,'\');
            save(fileName,'candidate');
                        
        end
        

    end
end

