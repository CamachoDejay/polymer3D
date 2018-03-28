function [tHold] = getTh(IM,titleStr)
%GETBGSIG Calculation of tresholds using the help of the user
%   Detailed explanation goes here

% create figure and ask for BG and signal areas
figure(1)
clf
imagesc(IM)
axis image
title(titleStr)

drawnow;
waitfor(helpdlg('Choose background'));
rect = getrect;
rect = round(rect);
xi = rect(2);
xf = xi + rect(4);
yi = rect(1);
yf = yi + rect(3);
bg = IM(xi:xf,yi:yf);

waitfor(helpdlg('Choose signal'));
rect = getrect;
rect = round(rect);
xi = rect(2);
xf = xi + rect(4);
yi = rect(1);
yf = yi + rect(3);
sig = IM(xi:xf,yi:yf);
bgV = bg(:);
sigV = sig(:);
mBg = mean(bg(:));
mSig = mean(sig(:));

%%
[tHold] = Misc.tholdSigBg(bgV,sigV);
% old method based on the mean of the BG and signal values
% tHold = mBg+(mBg+mSig)/2;
end

