
clear all
close all
clc

fPath = '/Users/rafa/Documents/MATLAB/data/Marisa';

post = 'DigestON_Post-Ex_1.tif';

pre = 'DigestON_Pre-Ex_1.tif';

%%
preFName = [fPath filesep pre];
postFName = [fPath filesep post];
[preIm, postIm] = expMic.getPrePost(preFName, postFName);

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

%%
IM = preIm;

[tHold] = expMic.getTh(IM,'PRE expansion');

[PreC, PreR] = expMic.segCellSquare(IM,tHold, 'PRE expansion');

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


IM = postIm;
titleStr = 'POST expansion';
[tHold] = expMic.getTh(IM,titleStr);

[PostC, PostR] = expMic.segCellSquare(IM,tHold, titleStr);

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

%%
[expF1, regCell] = expMic.getExpFactor(PreC,PostC);
fprintf('Cell based linear expansion factor: %2.4f \n', expF1)

figure(1)
clf
subplot(1,2,1)
plot(PreC.contNorm(1,:),PreC.contNorm(2,:))
hold on
plot(regCell.contour(1,:),regCell.contour(2,:))
plot(PreR.contNorm(1,:),PreR.contNorm(2,:))
[tmp] = expMic.appReg(PostR.contNorm,regCell.rot,regCell.scaling);
plot(tmp(1,:),tmp(2,:))
hold off
axis image
title('Cell based image registration')

[expF2, regR] = expMic.getExpFactor(PreR,PostR);
fprintf('Squared based linear expansion factor: %2.4f \n', expF2)

subplot(1,2,2)
plot(PreR.contNorm(1,:),PreR.contNorm(2,:))
hold on
plot(regR.contour(1,:),regR.contour(2,:))
plot(PreC.contNorm(1,:),PreC.contNorm(2,:))
[tmp] = expMic.appReg(PostC.contNorm,regR.rot,regR.scaling);
plot(tmp(1,:),tmp(2,:))
hold off
axis image
title('Squared based image registration')
%%
return

%% NOT WORKING
return
fix = BWpre;

F = griddedInterpolant(double(fix));
[sx,sy] = size(fix);
expF = 3;
xq = (0:1/expF:sx)';
yq = (0:1/expF:sy)';
vq = uint8(F({xq,yq}));

figure(4)
subplot(1,2,1)
imagesc(fix)
axis image
subplot(1,2,2)
imagesc(vq)
axis image
title('Higher Resolution')


mov = vq;
fix = BWpost;

tformEstimate = imregcorr(mov,fix);

Rfixed = imref2d(size(fix));
movReg = imwarp(mov,tformEstimate,'OutputView',Rfixed);

figure(5)
subplot(1,3,1)
imagesc(fix)
axis image
subplot(1,3,2)
imagesc(mov)
axis image
subplot(1,3,3)
imagesc(movReg)
axis image

%%
[optim, metric] = imregconfig('Multimodal');
regi = imregister(mov,fix,'rigid',optim,metric); 

figure(6)
subplot(1,3,1)
imagesc(fix)
axis image
subplot(1,3,2)
imagesc(mov)
axis image
subplot(1,3,3)
imagesc(regi)
axis image

%% NOT WORKING
return
mov = uint8(BWpre);
fix = uint8(BWpost);

[optim, metric] = imregconfig('Multimodal');
regi = imregister(mov,fix,'Similarity',optim,metric); 
% figure(3)
% imshowpair(regi,fix,'montage')

figure(3)
subplot(1,3,1)
imagesc(fix)
axis image
subplot(1,3,2)
imagesc(mov)
axis image
subplot(1,3,3)
imagesc(regi)
axis image
