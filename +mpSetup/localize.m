function [consLoc,totLoc] = localize(vol, rTh, GLRTprops)
%LOCALIZE detects molecules in volumetric multi-plane data
%   vol: volumetric data to analyze coming form multiplane setup
%   rTH: threshold for consolidation of positions in pixels

assert(size(vol,3)==8,'multiplane volumetric data mustconsist of 8 planes')
assert(isstruct(GLRTprops),'properties for GLRT must be a structure')

delta   = GLRTprops.delta;
pxSize  = GLRTprops.pxSnm;
FWHMnm  = GLRTprops.FWHMnm;
chi2    = GLRTprops.chi2;
FWHMpix = FWHMnm / pxSize;

% generate a list of all localized molecules per im plane, together with a
% common list for all im planes
totLoc{8,1}=[];
consLoc = [];

for i = 1:8
    % get one imaging plane
    im = vol(:,:,i);
    % change to double is needed for GLRT
    im = double(im);
    % find molecules in the image
    [ pos, meanFAR ] = Localization.smDetection( im, delta, FWHMpix, chi2 );
    % store positions
    totLoc{i} = [pos,meanFAR];
    
    % now we consolidate the positions
    
    if isempty(consLoc)
        consLoc = pos;
    else
        if ~isempty (pos)
            [ consLoc ] = Localization.consolidatePos( consLoc, pos, rTh );
        end
        
    end
    
end
end

