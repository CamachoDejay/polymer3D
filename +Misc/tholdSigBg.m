function [tHold] = tholdSigBg(bg,sig)
%THOLDSIGBG finds optimal treshhold that separates two distributions one
%comming from background and another from signal
%   Detailed explanation goes here

sigV = sig(:);
bgV  = bg(:);
assert(mean(bgV)<mean(sigV),'strange your background looks larger than your signal')

% find CDF of signal and CCDF of background
[sigCDF] = Plotting.getCDF(sigV);
[~, bgCCDF]  = Plotting.getCDF(bgV);

% find optimal tHold that separates two distributions
allX = [sigCDF.x; bgCCDF.x];
minX = min(allX);
maxX = max(allX);
globXq = linspace(double(minX),double(maxX),length(allX)*3);
sigVq = interp1(double(sigCDF.x),double(sigCDF.y),globXq);
bgVq  = interp1(double(bgCCDF.x),double(bgCCDF.y),globXq);
diffVq = abs(sigVq-bgVq);
[~,idx] = nanmin(diffVq);
tHold = globXq(idx);

end

