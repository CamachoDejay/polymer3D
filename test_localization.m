% test of localization algorithm
clear
close all
clc


% Information about detector
x_size = 100;
pix_size = 105; %[nm/pix]

% information about emitters
em_n = 20;
em_mean_int = 10000;
em_int_sigma = 100;

% information about psf
FWHM_nm  = 350; %[nm]

% information about normal bg
mean_bg = 1000;
SNR     = 3;
% bg_FWHM = 100;
d_range   = 'uint16';

% calculations for emitters positions and int
[ em_pos ] = EmitterSim.getRandPos( x_size, em_n );

int_model.name     = 'normal';
int_model.mean_int = em_mean_int;
int_model.sigma    = em_int_sigma;
[ em_int ] = EmitterSim.getIntensity( int_model, em_n );

% small calculations for image generation
y_size = x_size;
xv    = 1:1:x_size;
yv    = 1:1:y_size;
[X,Y] = meshgrid(xv,yv);
im = uint16(zeros(size(X)));

% calculations for psf
FWHM_pix = FWHM_nm / pix_size; %[pix]
sigma_pix = FWHM_pix / (2*((2*log(2))^0.5));
psf_model.name = 'gaussian';
psf_model.sigma_x = sigma_pix;
psf_model.sigma_y = psf_model.sigma_x;

[ G ] = EmitterSim.getPSF( X, Y, round(mean(xv)), round(mean(yv)), psf_model);
G = G.*(em_mean_int);
max_int = round(max(G(:)));

bg_FWHM = round(max_int/SNR);

% calculations for normal bg
[ bg_im ] = ImageNoise.gauss_noise( size(X), mean_bg, bg_FWHM, d_range );


for em_i = 1:em_n
    % get possition and intensity of the emitter
    x0_pix = em_pos(em_i,1);
    y0_pix = em_pos(em_i,2);
    n_counts = em_int(em_i);

    % get psf pdf of emitter
    [ psf ] = EmitterSim.getPSF( X, Y, x0_pix, y0_pix, psf_model );
    % get out only area of interest
    [ roi_lims ] = EmitterSim.getROI(x0_pix, y0_pix, 20, x_size, y_size);
    psf_roi = psf(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2));
    % sample psf
    [ im_roi ] = EmitterSim.samplePSF( psf_roi, n_counts, false );
    % update total image
    im(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2)) = ...
              im(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2)) + im_roi;

end

% add poisson noise
ccd_frame = imnoise((im + bg_im),'poisson');

figure(1)
imagesc( ccd_frame )
axis image


%%
im_in = double(ccd_frame);
delta = 4;
% calculations for psf
FWHM_nm  = 350; %[nm]
pix_size = 105; %[nm/pix]
FWHM_pix = FWHM_nm / pix_size; %[pix]
% for GLRT
chi2 = 24;

[ pos, inten ] = Localization.smDetection(im_in, delta, FWHM_pix, chi2 );


figure(1)
imagesc(im_in)
hold on
scatter(pos(:,1),pos(:,2),'or')
hold off
axis image




