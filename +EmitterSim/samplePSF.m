function [ im ] = samplePSF( psf, n_counts, do_plots )
%SAMPLEPSF This function generates the image of an emitter with a given psf
%and number of counts detected.
%   psf: point spread funtion of the emitter.
%   n_counts: total number of counts detected for the emitter.
%   do_plot: if true the function generate figures.
%   I wish to pick positions in space with a probability given by the PSF.
%   A way to do this is to use tha "datasample" function of matlab with
%   weights. When weigths are used, "datasample" samples the given vector
%   with probability proportional to the elements of weights. Weights must
%   be a vector with the same number of elemets as data elements in the
%   input vector.

% first I make sure that the number of counts is an integer.
n_counts = uint32(n_counts);

% now note that the probability given by psf comes as a 2D matrix. Thus to
% use "datasample" I must change psf into a vector. At the same time I need
% to keep track of the [x,y] coordinate value for each element of this
% vector. That is done below. 

% vector containing unique identifiers, this will be used to later recover
% the row and column value.
indx_vec = uint32( 1:(size(psf,1)*size(psf,2)));
% generate vectors to keep track of the row and colom value of each
% identifier
row_mat = uint32((1:size(psf,1))'* ones(1,size(psf,2)));
col_mat = uint32(((1:size(psf,2))'* ones(1,size(psf,1)))');
row_vec = row_mat(:);
col_vec = col_mat(:);
% generate weights for data sample.
w       = psf(:);
% sample the pixels of the iamge with a weight given by the psf
y       = datasample(indx_vec,n_counts,'Weights',w);
% find coordinates for each count
row = row_vec(y);
col = col_vec(y);
% init output image
im = uint16(zeros(size(psf)));
% populate the image with counts

% % im2 = uint32(zeros(size(psf)));
% % This approach looks smart but is very slow, is faster to do the for loop
% % tic
% % coors = [row(:) col(:)];
% % unique_pixels = unique(coors,'rows');
% % for ii=1:length(unique_pixels)
% %     c_val = sum(ismember(coors,unique_pixels(ii,:),'rows'));
% %     im2(unique_pixels(ii,1),unique_pixels(ii,2)) = c_val;
% % end
% % toc

 for i = 1:n_counts
    im(row(i),col(i)) = im(row(i),col(i))+1;
 end

 % if the user is interested now we create output images.
if do_plots
    figure(11)
    set(gcf,'position',[1290 140 380 780])
    subplot(2,1,1)
    imagesc(psf)
    axis image
    colorbar
    colormap 'gray'
    colormap(flipud(colormap))
    title('Initial PSF')

    subplot(2,1,2)
    imagesc(im)
    axis image
    colorbar
    colormap 'gray'
    colormap(flipud(colormap))
    title(['Emitter image for ' num2str(n_counts,'%1.1e') ' counts'])

    figure(12)
    set(gcf,'position',[740 80 450 840])
    subplot(2,1,1)
    surf(psf,'EdgeColor','none') 
    xlim([1 size(psf,1)])
    ylim([1 size(psf,2)])
    title('Initial Smooth PSF')
    
    subplot(2,1,2)
    surf(im,'EdgeColor','none') 
    xlim([1 size(im,1)])
    ylim([1 size(im,2)])
    title(['Detection events for ' num2str(n_counts,'%1.1e') ' counts'])
    
end

end

