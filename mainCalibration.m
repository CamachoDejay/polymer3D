% this script is a work in progress, the plan is to fix here the quick
% callibration of the channels. Basically I want to be able to lead the
% movies and built a dataset of [xDim,yDim,channel,frame] where all
% channels represent the same region of the image to a pixel resolution.
clear
close all
clc

fPath = '/Users/rafa/Documents/MATLAB/data/multi-plane/Beads - MC_4';
fName = 'Beads - MC_4_MMStack_Pos0.ome.tif';

fPath = [fPath filesep fName];

[frameInfo, movInfo, tnf ]= loadMovie.ome.getInfo( fPath );

posInfo = cat(1,frameInfo.Pos);

[ movC1, movC2, idx ] = loadMovie.ome.load( frameInfo, movInfo, 1:movInfo.maxFrame );

cal.flipc2 = true;
if cal.flipc2
    movC2 = flip(movC2,2);
end
%%


% ave_im1 = mean(movC1,3);
% ave_im2 = mean(movC2,3);

max_im1 = max(movC1,[],3);
max_im2 = max(movC2,[],3);

[ ch_c1, ~, common_w1 ] = mpSetup.cali.findChannels( max_im1 );
[ ch_c2, ~, common_w2 ] = mpSetup.cali.findChannels( max_im2 );

c_win = min([common_w1; common_w2]);
im_size = size(max_im1);
[ cal.ROI ] = mpSetup.cali.defineROI( c_win, ch_c1, ch_c2, im_size );

figure(1)
subplot(2,1,1)
imagesc(max_im1)
axis image
for i = 1:4
    rectangle('Position', cal.ROI(i,:))
end
subplot(2,1,2)
imagesc(max_im2)
axis image
for i = 5:8
    rectangle('Position', cal.ROI(i,:))
end

% chData has dim im_size1 im_size2 4channels Nframes
[ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI );

Z1 = posInfo(idx(:,1),3);
Z2 = posInfo(idx(:,2),3);
[ focus_met, in_focus ] = mpSetup.cali.getFocusMetric( chData1c, chData2c , Z1, Z2 );
figure(2)
for i = 1:4
    plot(Z1,focus_met(:,i))
    hold on
end
for i = 5:8
    plot(Z2,focus_met(:,i))
end
hold off

[ cal.neworder, in_focus ] = mpSetup.cali.getNewOrder( in_focus );


data12 = cat(3,chData1c,chData2c);

[ im_shifts ] = mpSetup.cali.simpleImShift( in_focus, data12 );

[ cal.ROI ] = mpSetup.cali.refineROI( cal.ROI, im_shifts );

% chData has dim im_size1 im_size2 4channels Nframes
[ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI );

data12 = cat(3,chData1c,chData2c);

[ ~, cal.Icorrf ] = mpSetup.cali.findChInt( data12, in_focus );


[ data12 ] = mpSetup.cali.arrange_correctInt( data12, cal );

return
%% just to show that it works

df = 20;
pos = [130,170];
imagesc(data12(pos(2)-df:pos(2)+df,pos(1)-df:pos(1)+df,1,25))
axis image

testData = data12(pos(2)-df:pos(2)+df,pos(1)-df:pos(1)+df,:,40);

for i = 1:8
    subplot(2,4,i)
    imagesc(testData(:,:,i))
    axis image
end
%% TODO I still have to fix this part, however I do not think is needed for our purposes

% init output
    tf=cell(sys.nplanes-1,1);

    for m=(1:sys.nplanes-1)-1
        
        ch = sys.nplanes-m

        % get image sequence for 2 channels to compare
        img_seq_ch_fix=squeeze(data12(:,:,ch,:));
        img_seq_ch_moving=squeeze(data12(:,:,ch-1,:));

        % see when each channel is in focus % TODO this I can do better
        focus_metric_fix = squeeze(mean(max(img_seq_ch_fix)));
        focus_metric_moving = squeeze(mean(max(img_seq_ch_moving)));
        
        % get a frame when both cameras are simultaneously in best focus
        [~,indcom]=max(focus_metric_fix.*focus_metric_moving);
        im_fix=medfilt2(img_seq_ch_fix(:,:,indcom),[3 3]);
        im_mov=medfilt2(img_seq_ch_moving(:,:,indcom),[3 3]);
        
        % identify corresponding beads in the two channels
        [my, mx] = ccrShiftEstimation(im_fix,im_mov,10); % ! Stefan switched mx and my
        

        % extracting the center of gravities of the beads
        out=struct;
        out.ru=5;
        sys.bg=1;
        out.nh=0;

        out = hriSegmentation(double(im_fix),cal.bgth,cal.logsize,out);
        out = hriFilterSegments(double(im_fix),cal.aupl,cal.alol,sys,out);
        cog_fix = [out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];
        
        out = hriSegmentation(double(im_mov),cal.bgth,cal.logsize,out);
        out = hriFilterSegments(double(im_mov),cal.aupl,cal.alol,sys,out);
        cog_mov = [out.xoAbsolute+0.1*out.xoEstimate, out.yoAbsolute+0.1*out.yoEstimate];

        %identify corresponding center of gravities
        
        pixel_tolerance = 2;
        mainSet = cog_fix - repmat([mx my],size(cog_fix,1),1);
        testSet = cog_mov;
        
        [ cog_common_fix, cog_common_mov ] = ...
             ImageProc.consolidatePos( mainSet, testSet, pixel_tolerance );
         
        d = cog_common_mov - cog_common_fix;
        
        disp(['avg error of COG coordinates using pure displacement at pixel level: ' num2str(mean(sqrt(d(:,1).^2+d(:,2).^2)))]);
        
        cog_common_fix     = cog_common_fix + ...
                                 repmat([mx my],size(cog_common_fix,1),1);
        
        
        % affine transformation
        ccm=fliplr(cog_common_mov);
        ccf=fliplr(cog_common_fix);
        % Infer spatial transformation from control point pairs
        tf{m+1}=cp2tform(ccm,ccf,'similarity');
        
%         % I think this is here for example purposes
%         % Apply forward spatial transformation.
%         [x,y] = tformfwd(tf{m+1},cog_common_mov(:,2),cog_common_mov(:,1));
    end
    
    % tf{1} transforms ch7 into ch8
    % tf{2} transforms ch6 into ch7
    % tf{3} transforms ch5 into ch6
    % tf{4} transforms ch4 into ch5
    % tf{5} transforms ch3 into ch4
    % tf{6} transforms ch2 into ch3
    % tf{7} transforms ch1 into ch2
       
    cal.tf = tf;

