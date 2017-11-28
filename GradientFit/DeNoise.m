% perform de-noising by using band pass filter
% -------------------------------------------------------------------------
% Inputs
%   im    : Raw image
% -------------------------------------------------------------------------
% Outputs
%   imDN  : Denoised image
% -------------------------------------------------------------------------

function imDN = DeNoise(im)

[X,Y] = meshgrid(-4:4);
S1 = 0.5;
S2 = 1.8;

G1 = exp(-(X.^2+Y.^2)/(2*pi*S1^2));
G2 = abs(1-exp(-(X.^2+Y.^2)/(2*pi*S2^2)));

BandPass = G1/sum(G1(:)) - G2/sum(G2(:));

imDN = conv2(im,BandPass,'same');

end