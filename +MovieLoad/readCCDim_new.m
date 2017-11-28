function [imFrame, datalength, info] = readCCDim_new(filename,infor,fr)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CopyRight Hiroshi Ujii & Beno?t Muls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function load fr-th frame of spe file
% filename is the location of the spe file you want to load
% fr is number of frame you want to read
% if infor=1 to know data length
% if infor=2 to get the matrix
% type=1:SPE
% type=2:sif
% type=3: tiff, tif
% type=4; bmp, gif, jpeg
% type=5; his
%
%Output
% imFrame; image (when infor==2)
% datalength; number of frame (when infor==3)
% info.extme
% info.data
    %% 
    if nargin<2
        infor=1;
        fr=1;
    elseif nargin<3
        fr=1;
    end

    %%
    [~,~,ext]=fileparts(filename);
    switch lower(ext)
        case '.spe'  % SPE
            type=1;
            skp=4100;
        case '.sif'  % Andor CCD camera
            type=2;
        case {'.tif' '.tiff'}
            type=3;
        case {'.bmp' '.gif' '.jpg' '.jpeg' '.png'}
            type=4;
        case '.his'  % Hamamatsu CCD camera
            type=5;
        otherwise %unknown format
            imFrame=[];datalength=[]; info=[];
            return;
    end
    %%
    info.extme='not available';
    info.date='not available';
    switch type
        case 1 %SPE
            spefile = fopen(filename,'r','ieee-le');
            fseek(spefile, 10, -1);
            [extme, count]=fread(spefile, 1, 'float32');%exposure time in second
            extme=extme*1000; % exposure time in ms
            info.extme=[num2str(extme) 'ms'];
            fseek(spefile, 20, -1);
            [mdate, count]=fread(spefile, 10, 'int8');
            info.date=char(mdate'); %measurement date
            fseek(spefile, 42, -1);
            [xsize, count]=fread(spefile, 1, 'uint16');%actual # of pixels on x axis
            fseek(spefile, 656, -1);
            [ysize, count]=fread(spefile, 1, 'uint16');%actual # of pixels on y axis
            fseek(spefile, 1446, -1);
            [datalength, count]=fread(spefile, 1, 'uint16');%number of frames in file
            fseek(spefile, 108, -1);
            [datatype, count]=fread(spefile, 1, 'int16');%datatype 3*unit16, 2:int16, 1:int32, 0:float32
            fseek(spefile, 150, -1);
            [BackGrndApplied, count]=fread(spefile, 1, 'uint16'); %1: background subdone
            fclose(spefile);
        case 2 %Andor .sif file
            fid=fopen(filename, 'r');
            fseek(fid, 0,1);
            info.filesize=ftell(fid);
            frewind(fid);
            fgetl(fid);
            fgetl(fid);
            line1=fgetl(fid);
            if length(line1)<100
                line1=[line1, ' ', fgetl(fid)];
            end
            infos=regexp(line1, '(\S)+\s', 'match');
            info.timeStamp=sscanf(infos{5}, '%f');
            info.timeStamp=datenum([1970 1 1 0 0 info.timeStamp]);
            info.time=datestr(info.timeStamp);
            info.temp=sscanf(infos{48}, '%f');
            info.exposureTime=sscanf(infos{13}, '%f');
            info.cycleTime=sscanf(infos{15}, '%f');
            info.EM=sscanf(infos{32}, '%f');
            version1 = sscanf(infos{55}, '%f');
            version2 = sscanf(infos{56}, '%f');
            version3 = sscanf(infos{57}, '%f');
            version4 = sscanf(infos{58}, '%f');
            info.version=sprintf('%d.%d.%d.%d', version1, version2, version3, version4);

            line2=fgetl(fid);
            detectors=regexp(line2, '(\S)+\s', 'match');
            info.detectorType=detectors{1};
            if strcmp(info.detectorType, 'DU897_BV')
                info.currentBitDepth = 14;
            else
                info.currentBitDepth = 11;
            end
            fgetl(fid);
            info.filename=fgetl(fid);

            fgetl(fid); fgetl(fid);
            headsize=2400;
            s=fread(fid, headsize, 'uint8=>char')';
            pos=strfind(s, 'Pixel number');
            fseek(fid, pos(end)-1-headsize, 0);
            s=fgetl(fid);
            s=regexp(s, '(\d)+\S', 'match');

            info.pixelsinframe=sscanf(s{end},'%d');
            info.length=floor(info.filesize/info.pixelsinframe/4);
            datalength=info.length;

            s=regexp(fgetl(fid), '(\S)+\s', 'match');
            info.region=[sscanf(s{2},'%d'), sscanf(s{4},'%d'), sscanf(s{5},'%d'), sscanf(s{3},'%d')];
            info.bin=[sscanf(s{6},'%d'), sscanf(s{7},'%d')];
            info.width=(info.region(2)-info.region(1)+1)/info.bin(1);
            info.height=(info.region(4)-info.region(3)+1)/info.bin(2);
            
            if strcmp(info.version(1:4), '4.21')
                skp=info.filesize-info.pixelsinframe*info.length*4-6;
            elseif strcmp(info.version(1:4), '4.24') && info.length==1
                skp=info.filesize-info.pixelsinframe*info.length*4-6-412;
            elseif strcmp(info.version(1:4), '4.24') && info.length>1
                skp=info.filesize-info.pixelsinframe*info.length*4-6-274;
            else
                skp=info.filesize-info.pixelsinframe*info.length*4-6;
            end
            fclose(fid);
        case 3 
            if  infor==1
                datalength=length(imfinfo(filename));
                %info.extme=[];
                % info.data=[];
            elseif infor==2
                datalength=1;
                %  info.extme=[];
                %info.data=[];
            end
        case 4
            datalength=1;
            %   info.extme=[];
            %info.data=[];
            fr=1;
        case 5 % Hamamatsu
            spefile = fopen(filename,'r','ieee-le');
            [header,count]=fread(spefile,32,'uint16');
            datalength=header(8);
            ysize=header(3);
            xsize=header(4);
            F=fread(spefile, 800, 'uint8');
            F=char(F');
            a=findstr(F, '~Hokawo~');
            if isempty(a);a=findstr(F, '~WASABI~');end
            skp=a+7;
            extme=findstr(F, 'vExpTim1');
            extme=str2num(F(extme+9:extme+16))*1000;% integration time in ms
            info.extme=[num2str(extme) 'ms'];
            date=findstr(F, 'vDate');
            date=F(date+8:date+15);
            info.date=date;
            vTStamp=findstr(F, 'vTStamp');
            vTStamp=str2num(F(vTStamp+8:vTStamp+22));% Tstamp
            info.TStamp=vTStamp;
            if strfind(F, 'vrcAcqArea')
                version=2;
            end
            fclose(spefile);
    end

    %%
    if infor == 1; imFrame=info; return; end

    if fr > datalength && type~=3
        errordlg('Frame number bigger than data length...')
        imFrame=[];
        return;
    end
    %% infor==2
    switch type
        case 1 
            %% .spe file
            spefile = fopen(filename,'r','ieee-le');
            if datatype==3;
                if BackGrndApplied==1
                    st=fseek(spefile, xsize*ysize*(fr-1)*4+skp,-1);
                else
                    st=fseek(spefile, xsize*ysize*(fr-1)*2+skp,-1);
                end
                %st=fseek(spefile, xsize*ysize*(fr-1)*2+skp,-1);
                if st==0;
                    imFrame=fread(spefile, xsize*ysize, 'uint16');
                    imFrame=double(imFrame);
                else
                    disp('error!!!')
                    return
                end
            elseif datatype ==2;
                if BackGrndApplied==1
                    st=fseek(spefile, xsize*ysize*(fr-1)*4+skp,-1);
                else
                    st=fseek(spefile, xsize*ysize*(fr-1)*2+skp,-1);
                end
                %st=fseek(spefile, xsize*ysize*(fr-1)*2+skp,-1);
                if st==0;
                    imFrame=fread(spefile, xsize*ysize, 'int16');
                    imFrame=double(imFrame);
                else
                    disp('error!!!')
                    return
                end
            elseif datatype ==1
                %st=fseek(spefile, xsize*ysize*(fr-1)*2+skp,-1);
                if BackGrndApplied==1
                    st=fseek(spefile, xsize*ysize*(fr-1)*4+skp,-1);
                else
                    st=fseek(spefile, xsize*ysize*(fr-1)*2+skp,-1);
                end
                if st==0;
                    imFrame=fread(spefile, xsize*ysize, 'int32');
                    imFrame=double(imFrame);
                else
                    disp('error!!!')
                    return
                end
            elseif datatype==0
                if BackGrndApplied==1
                    st=fseek(spefile, xsize*ysize*(fr-1)*4+skp,-1);
                else
                    st=fseek(spefile, xsize*ysize*(fr-1)*2+skp,-1);
                end
                if st==0;
                    imFrame=fread(spefile, xsize*ysize, 'float32');
                    imFrame=double(imFrame);
                else
                    disp('error!!!')
                    return
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%remove header and restructure data as 2D matrix%%%%%%%
            imFrame=reshape(imFrame, xsize,ysize);
            imFrame=imFrame';
            fclose(spefile);
        case 2 %andor .sif file
            spefile = fopen(filename,'r');   %'ieee-le'
            fseek(spefile, skp+info.width*info.height*(fr-1)*4,-1);
            imFrame=fread(spefile, [info.width, info.height],'uint16', 2);
            imFrame=flipud(imFrame');
            imFrame=double(imFrame);
            fclose(spefile);
            %imagesc(imFrame'); colormap hot
        case 3 %tif
            %% 
            imFrame=imread(filename, fr);
            if length(size(imFrame))==3
                imFrame=rgb2gray(imFrame);
            end
            % datalength=length(imfinfo(filename));
        case 4 %jpg,png,...
            %% 
            imFrame=imread(filename);
            if length(size(imFrame))==3
                imFrame=rgb2gray(imFrame);
            end
        case 5%Hokao
            %% 
            spefile = fopen(filename,'r','ieee-le');   %'ieee-le'
            if exist('version', 'var') && version==2
                str2find = '~Hokawo~';
                frewind(spefile);
                c = fgetl(spefile);
                b = strfind(c, str2find);
                
                if ~isempty(b)
                    head_len = b + length(str2find)-1;                    
                else
                    frewind(spefile);
                    head_len = length(str2find);
                    cc = fread(spefile, head_len, 'uint8=>char')';
                    r = true;
                    start_idx = [];
                    end_idx   = [];
                    max_hl    = 2000;
                    while r
                       cc = fread(spefile,  1, 'uint8=>char');
                       head_len = head_len + 1;
                       if strcmp(cc,'~')                       
                           if isempty (start_idx)
                                start_idx = head_len;
                           elseif isempty (end_idx)
                               end_idx    = head_len;
                               assert(end_idx-start_idx == 7,'Problems with header')
                               r = false;
                           end                       
                       else
                           r = head_len < max_hl;
                       end
                    end
                    assert(~isempty(end_idx),'Problems with header')
                end
%                 assert(~isempty(b),'Problems with header in HIS')
                frewind(spefile);

                for i=1:fr
                    c=fread(spefile,  head_len, 'uint8=>char')';
                    b=strfind(c, str2find);
                    assert(~isempty(b),'Problems with header in HIS')
                    if mod(b,2)==1
                        fseek(spefile,  b-head_len+7, 0);
                    else
                        fseek(spefile,  b-head_len+8, 0);
                    end
                    if i~=fr
                        fseek(spefile,ysize*xsize*2, 0);
                    else
                        imFrame=fread(spefile,[ysize,xsize], 'uint16');
                    end
                end
            else
                sk=32*fr*2 + skp*(fr)*1 + (xsize*ysize)*2*(fr-1);
                st=fseek(spefile, sk,-1);
                %st=fseek(spefile, skp+(xsize*ysize*(fr-1)*4),-1);
                imFrame=fread(spefile, xsize*ysize,'uint16');
                imFrame=reshape(imFrame, ysize,xsize);
                %imFrame=flipud(imFrame');
            end
            imFrame=imFrame';
            imFrame=double(imFrame);
            if 1
                [header,count]=fread(spefile,32,'uint16');
                F=fread(spefile, 500, 'uint8');
                F=char(F');

                vTStamp=findstr(F, 'vTStamp');
                vTStamp=str2num(F(vTStamp+8:vTStamp+22));% Tstamp in ms
                info.TStamp=vTStamp;
            end
            fclose(spefile);
    end
    imFrame=double(imFrame);
end
%%
function info=readinfo(filename)
    fileinfo=dir(which(filename));
    info.filesize=fileinfo.bytes;
    
    [~,~,ext]=fileparts(filename);
    switch lower(ext)
        case '.sif' % Andor v4.21.30006.0
            fid=fopen(filename, 'r');
            fgetl(fid); fgetl(fid);
            line1=fgetl(fid);
            infos=regexp(line1, '(\S)+\s', 'match');
            info.timeStamp=sscanf(infos{5}, '%f');
            info.timeStamp=datenum([1970 1 1 0 0 info.timeStamp]);
            info.time=datestr(info.timeStamp);
            info.temperature=sscanf(infos{48}, '%f');
            info.exposureTime=sscanf(infos{13}, '%f');
            info.cycleTime=sscanf(infos{15}, '%f');
            info.EM=sscanf(infos{32}, '%f');
            version1 = sscanf(infos{55}, '%f');
            version2 = sscanf(infos{56}, '%f');
            version3 = sscanf(infos{57}, '%f');
            version4 = sscanf(infos{58}, '%f');
            info.version=sprintf('%d.%d.%d.%d', version1, version2, version3, version4);

            line2=fgetl(fid);
            detectors=regexp(line2, '(\S)+\s', 'match');
            info.detectorType=detectors{1};
            if strcmp(info.detectorType, 'DU897_BV')
                info.currentBitDepth = 14;
            else
                info.currentBitDepth = 11;
            end
            fgetl(fid);
            info.filename=fgetl(fid);

            fgetl(fid); fgetl(fid);
            headsize=2400;
            s=fread(fid, headsize, 'uint8=>char')';
            pos=strfind(s, 'Pixel number');
            fseek(fid, pos(end)-1-headsize, 0);
            s=fgetl(fid);
            s=regexp(s, '(\d)+\S', 'match');

            if length(s)==5
                info.length=1;
                info.pixelsinframe=sscanf(s{5},'%d');
            elseif length(s)==6
                info.length=sscanf(s{4},'%d');
                info.pixelsinframe=sscanf(s{6},'%d');
            end

            s=regexp(fgetl(fid), '(\S)+\s', 'match');
            info.region=[sscanf(s{2},'%d'), sscanf(s{4},'%d'); sscanf(s{5},'%d'), sscanf(s{3},'%d')];
            info.bin=[sscanf(s{6},'%d'), sscanf(s{7},'%d')];
            info.width=(info.region(1,2)-info.region(1,1)+1)/info.bin(1);
            info.height=(info.region(2,2)-info.region(2,1)+1)/info.bin(2);

            info.skip=info.filesize-info.pixelsinframe*info.length*4-6;
            fclose(fid);
        case {'.his', '.rbf'} % 
            fid=fopen(filename, 'r');
            header0=fread(fid,32,'uint16');
            info.version=header0(19);

            header=fgetl(fid);
            vSens=regexp(header, 'vSens=([^;]*);', 'match', 'once');
            info.EM=str2double(vSens(7:end-1));
            vDate=regexp(header, 'vDate=([^;]*);', 'match', 'once');
            info.date=vDate(7:end-1);
            vTStamp=regexp(header, 'vTStamp=([^;]*);', 'match', 'once');
            info.timeStamp=str2double(vTStamp(9:end-1));
            vExpTim1=regexp(header, 'vExpTim1=([^;]*);', 'match', 'once');
            info.exposureTime=str2double(vExpTim1(10:end-1));
            temperature=regexp(header, 'vCTemp=([^;]*);', 'match', 'once');
            info.temperature=str2double(temperature(8:end-1));

            vrcSub=regexp(header, 'vrcSub=([^;]*);', 'match', 'once');
            vrcSub=vrcSub(8:end-1);
            info.region=sscanf(vrcSub, '%d, %d, %d, %d');
            info.region=[info.region(1)+1,  info.region(3); info.region(2)+1,  info.region(4)];
            binX=regexp(header, 'vBinX=\d', 'match', 'once');
            info.bin(1)=str2double(binX(7:end));
            binY=regexp(header, 'vBinY=\d', 'match', 'once');
            info.bin(2)=str2double(binY(7:end));
            info.width=(info.region(1,2)-info.region(1,1)+1)/info.bin(1);
            info.height=(info.region(2,2)-info.region(2,1)+1)/info.bin(2);

            if strcmpi(ext, '.his')
                info.skip0=64+regexp(header, '~Hokawo~', 'end');
                frewind(fid);    fseek(fid, info.skip0+512*512*2, 0);
                headerFrame=fgetl(fid);
                info.skip=regexp(headerFrame, '~Hokawo~', 'end');
            elseif strcmpi(ext, '.rbf')
                info.skip=64+regexp(header, '~Hokawo~', 'end');
            end
            fclose(fid);
        case {'.bmp' '.jpg' '.jpeg' '.png'}
            info=imfinfo(filename);
            info.length=1;
        case {'.tif' '.tiff' '.gif'}
            infos=imfinfo(filename);
            info=infos(1);
            info.length=length(info);
        case {'.avi' '.mpg' '.wmv' '.mp4' '.mov'}
            info=VideoReader(filename);
            info.length=info.NumberOfFrames;
    end
end
