function [ totCor, Icor ] = findChInt( cam1, cam2, inFocus, doFigure )
%FINDCHINT finds hte intensity difference between the different channels
%   Detailed explanation goes here

switch nargin
    case 3
        doFigure = false;
end

% number of frames
N = size(cam1,4);
% intensity per channel
I = zeros(size(cam1,3)+size(cam2,3),N);
% calculating mean int of each channel for each frame
I(1:size(cam1,3),:) = mean(mean(cam1,2),1);
I(size(cam1,3)+1:size(cam1,3)+size(cam2,3),:) = mean(mean(cam2,2),1);
   
% to avoid problems with bleaching I can calculate Int differences between
% channels by using the frame where both channels are almost in focus. 
% global channel order
globCh = cat(1,inFocus.globalch);
% frame of sharp focus
focus  = cat(1,inFocus.frame);
% find the order of the list
[~, idx] = sort(globCh);
% get grames of sharp focus
focus = focus(idx);
% calculate the common frame to use
dF = diff(focus);
commonF = round(focus(1:end-1)+dF./2);

% calculate intensity correction
% Icor contains the reference channel, the channel to transform, the frame
% used and the correction factor.
Icor = ones(length(focus),4);
Icor(1,1) = idx(1);
Icor(1,2) = idx(1);
for i = 2:length(focus)
    % index of the first channel
    idx1 = idx(i-1);
    % index of the second channel
    idx2 = idx(i);
    % frame used 
    F    = commonF(i-1);
    % mean intensity values
    I1 = I(idx1,F);
    I2 = I(idx2,F);
    % now we store all values
    Icor(i,1) = idx1;
    Icor(i,2) = idx2;
    Icor(i,3) = F;
    % I must multiply by this factor to change int of channel idx2 into the
    % same level of idx1;
    Icor(i,4) = I1/I2;
end

% so far the corrections are between different channels, I have to find a
% factor to correct all to the same refence. I do this by its cumulative
% multiplication
totCor = cumprod(Icor(:,4));
% now I just have to order it correctly
[~,newOrd] = sort(Icor(:,2));
% sorting
totCor = totCor(newOrd);

if doFigure
    Itest = I.*repmat(totCor,1,N);


    figure(1)
    subplot(1,2,1)
    plot(1:N,I(1,:))
    hold on
    for i=2:8
       plot(1:N,I(i,:)) 
    end
    hold off
    xlim([1, N])

    subplot(1,2,2)
    plot(1:N,Itest(1,:))
    hold on
    for i=2:8
       plot(1:N,Itest(i,:)) 
    end
    hold off
    xlim([1,N])

end

% old method
% % normalizing the plots to its max value 
% Inorm = I./ repmat(max(I,[],2),1,N);
% % using normalized traces we can have an idea about the bleaching in the
% % sample
% bleaching = mean(Inorm,1);
% % correction for the blaching
% Icor = I./ repmat(bleaching,8,1);
%     
% Ich = zeros(8,1);
% for i = 1:8
%     Ich(i) = Icor(i,inFocus(i).frame);
% end
% 
% Ich_norm = max(Ich)./Ich;

end

