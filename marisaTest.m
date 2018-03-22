clear
close all
clc

%%%%% USER INPUT %%%%%%%%%%%%%
% path to folder contining the pre and post expanded images in tif format
fPath = '/Users/rafa/Documents/MATLAB/data/Marisa';

% name of the pre and post expanded images, .tif!!!
pre = 'DigestON_Pre-Ex_1.tif';
post = 'DigestON_Post-Ex_1.tif';

%%%%%% END OF USER INPUT %%%%%%%

%% load and crop ROIs
% full paths to files
%   pre
preFName = [fPath filesep pre];
%   post
postFName = [fPath filesep post];
%   get ROIs from the images
[preIm, postIm] = expMic.getPrePost(preFName, postFName);

% asking user if all goes ok
answer = questdlg('Should we continue?', ...
	'user test', ...
	'Yes','No','Yes');
% Handle response
switch answer
    case 'Yes'
        disp('we continue')
    case 'No'
        return
    otherwise
        error('ASDASD')
end

%% segment cell and square
% pre expansion
IM = preIm;
% get image threshold
titleStr = 'PRE expansion';
[tHold] = expMic.getTh(IM,titleStr);
% do segmentation
[PreC, PreR] = expMic.segCellSquare(IM,tHold, titleStr);
% ask is all goes well
answer = questdlg('Should we continue?', ...
	'user test', ...
	'Yes','No','Yes');
% Handle response
switch answer
    case 'Yes'
        disp('we continue')
    case 'No'
        return
    otherwise
        error('ASDASD')
end

% post expansion
IM = postIm;
% find threshold
titleStr = 'POST expansion';
[tHold] = expMic.getTh(IM,titleStr);
% segment cell and square
[PostC, PostR] = expMic.segCellSquare(IM,tHold, titleStr);
% ask if all goes well
answer = questdlg('Should we continue?', ...
	'user test', ...
	'Yes','No','Yes');
% Handle response
switch answer
    case 'Yes'
        disp('we continue')
    case 'No'
        return
    otherwise
        error('ASDASD')
end

%% calculating the linear exp factor based on the contours
% cell contour based registration and expansion factor
[expF1, regCell] = expMic.getExpFactor(PreC,PostC);
[expF2, regSq] = expMic.getExpFactor(PreR,PostR);

figure(1)
clf
% image of pre expanded contour - FIX
p1 = plot(PreC.contNorm(1,:),PreC.contNorm(2,:),'g','linewidth',3);
hold on
% image of post expanded registered image
p2 = plot(regCell.contour(1,:),regCell.contour(2,:), ':k','linewidth',2);
% image of pre expanded squared - FIX
plot(PreR.contNorm(1,:),PreR.contNorm(2,:),'g','linewidth',3)
% image of post expanded registered square 
[tmp] = expMic.appReg(PostR.contNorm,regCell.rot,regCell.scaling);
plot(tmp(1,:),tmp(2,:),':k','linewidth',2)
hold off
axis image
title('Cell based image registration')
a = gca;
a.FontSize = 14;
legend([p1 p2],'Pre expansion','Post registered',...
                      'Location','northoutside','Orientation','horizontal')

figure(2)
clf
% image of pre expanded square - FIX
p1 = plot(PreR.contNorm(1,:),PreR.contNorm(2,:),'g','linewidth',3);
hold on
% image of post expanded registered square
p2 = plot(regSq.contour(1,:),regSq.contour(2,:),':k','linewidth',2);
% image of pre expanded cell - FIX
plot(PreC.contNorm(1,:),PreC.contNorm(2,:),'g','linewidth',3)
[tmp] = expMic.appReg(PostC.contNorm,regSq.rot,regSq.scaling);
% image of post expanded registered cell
plot(tmp(1,:),tmp(2,:),':k','linewidth',2)
hold off
axis image
title('Squared based image registration')
a = gca;
a.FontSize = 14;
legend([p1 p2],'Pre expansion','Post registered',...
                      'Location','northoutside','Orientation','horizontal')

fprintf('-- output ---- \n')
fprintf('Cell based registration: \n')
fprintf('Cell based linear expansion factor: %2.4f \n', expF1)
fprintf('Cell based registration mean error: %2.4f [pre pix]\n', regCell.error)
fprintf('------ \n')
fprintf('Square based registration: \n')
fprintf('Square based linear expansion factor: %2.4f \n', expF2)
fprintf('Square based registration mean error: %2.4f [pre pix]\n', regSq.error)
fprintf('---- end of output ---- \n')

