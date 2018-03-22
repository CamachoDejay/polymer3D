function [preIm, postIm] = getPrePost(preFPath, postFPath)
%GETPREPOST small code to crop the cell of interest from pre and post
%expanded images
%   Detailed explanation goes here

% try to load, if not you have to give me a tif
try
   [ preIm ] = loadMovie.tif.getframes( preFPath, 1 );
catch ME
   error('Input for pre exp must be a tif')
end 

try
   [ postIm ] = loadMovie.tif.getframes( postFPath, 1 );
catch ME
   error('Input for pre exp must be a tif')
end 

% simple crop of the images
IM = preIm;
figure(1)
clf
imagesc(IM)
title('Pre expansion image')
axis image
drawnow;
waitfor(helpdlg('Choose cell to process'));
rect = getrect;
rect = round(rect);
xi = rect(2);
xf = xi + rect(4);
yi = rect(1);
yf = yi + rect(3);
preIm = IM(xi:xf,yi:yf);

IM = postIm;
figure(1)
clf
imagesc(IM)
title('Post expansion image')
axis image
drawnow;
waitfor(helpdlg('Choose cell to process'));
rect = getrect;
rect = round(rect);
xi = rect(2);
xf = xi + rect(4);
yi = rect(1);
yf = yi + rect(3);
postIm = IM(xi:xf,yi:yf);

% output figure
figure(1)
clf
subplot(1,2,1)
imagesc(preIm)
title('Pre expansion image to process')
axis image
drawnow;
waitfor(2)

subplot(1,2,2)
imagesc(postIm)
title('Post expansion image to process')
axis image
drawnow;
waitfor(2)
end

