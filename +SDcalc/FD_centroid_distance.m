function [ FD, FDnames ] = FD_centroid_distance( r, t, n_comp, do_figure )
% For a given shape defined by a closed curve C and tha is represented by a
% one dimentional function called shape signature we can define fourier
% descriptors. In this case we use centroid distance as shape signature.
%   To achieve translation, rotation invariance, phaseinformation of the
%   F_z(n) are ignored and only the magnitudes abs(F_z(n)) are used. Scale
%   invariance is achieved by dividing the magnitudes by the DC component,
%   i.e., abs(F_z(1)). The normalized Fourier coefficients are called FD.
%   The similarity between a query shape Q and a target shape T is measured
%   by the Euclidean distance between their FD representations.

% % raw_r     = centroid_distance.raw.r;
% % raw_t     = centroid_distance.raw.t;
% % raw_theta = centroid_distance.raw.theta;

FDnames = cell(n_comp,1);
for i = 1:n_comp
    FDnames(i) = {['FD_' num2str(i)]};        
end

if any([isempty(r),isempty(t)])
    % Generation of output
    FD.components = []; 
    FD.fourier    = [];

else

r     = r(:); %% IMPORTANT??? CHECK IF flipud is needed (rotation convention)
m     = length(r);          % Window length
% n     = pow2(nextpow2(m));  % Transform length
F_r   = fft(r,m);           % DFT
% power = F_r.*conj(F_r)/m;     % Power of the DFT
magnitudes = abs(F_r);           % Amplitude

if do_figure
    figure(3)
    subplot(2,1,1)
    plot(t,r,'b')
    title('cent dist func -r(t)')
    xlim([0 2*pi])

    F_zs = fftshift(F_r)./F_r(1);
    subplot(2,1,2)
    plot(real(F_zs),'b')
    hold on
    plot(imag(F_zs),'--g')
    hold off
    legend('real','imaginary')
    xlim([500-n_comp 502+n_comp])
     title('Fourier transform of r')
end
%%
% n_comp = 6; % I can use 60 comp to smoth the figure out and remove errors
             % on the other hand 4-6 comp seems to be enough to describe the
             % cells 
             
if do_figure             
    F_z2fit = (F_r);%./(F_z(1));
    F_z2fit (n_comp+2:end-n_comp) = 0;
%     FD = [F_r(end-n_comp+1:end); F_r(1:n_comp)];
%     FD_power = [power(end-n_comp+1:end); power(1:n_comp)];
%     FD_amp   = [ampli(end-n_comp+1:end); power(1:n_comp)];
    z_fit  = ifft (F_z2fit);

    figure(4)
    plot(r,'-k')
    hold on
    plot(z_fit,'--r')
    plot(r-z_fit,'-g')
    hold off
end
% F_z(1) is the DC component and its value scales with the size of the
% shape. Thus, I should normalize by it wen calculating the fourier
% descriptors.
magnitudes_norm = magnitudes./(magnitudes(1));
% Further, centroid distance is a real function, therefore, only half of
% the coefficients are distinct.

% Generation of output
FD.components = magnitudes_norm(1+1:n_comp+1); %remember DC component will be always 1
FD.fourier    = F_r;
FD.fourier(n_comp+2:end-n_comp) = 0;

end

% % % %     subplot(1,2,2)
% % % %     % if I use fit values this means little, as I dont know the angles anymore.
% % % %     % at this point in time I dont know of an easy way to recover the fit
% % % %     % angles theta, thus I can not reconstruct the images directly. However, I
% % % %     % have tested it with the raw data and I know it works.
% % % %     polar(raw_theta, raw_r ,'k')
% % % %     hold on
% % % %     polar(raw_theta, z_fit,'--r')
% % % %     hold off

end

