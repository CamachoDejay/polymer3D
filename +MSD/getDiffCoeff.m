%Function to calculate diffusion coefficent based on MSD curve

% MSD curve needs to be provided as a vector
% fitRange is a integer, how many data point to be consider for the fit
% dim is the dimension in which MSD was calculated (1D,2D or 3D)

function D = getDiffCoeff(msd,fitRange,dim)
    
    switch dim
        case '1D'
            div = 2;
        case '2D'
            div = 4;
        case '3D'
            div = 6;
        otherwise
            error('Unknown dim, dim needs to be provided as 1D 2D or 3D')
    end

    assert(min(size(msd))==1,'MSD needs to be provided as a vector')
    assert(and(fitRange<=1,isnumeric(fitRange)),'fit Range needs to be numerical between 0 and 1');
    
    tofit = msd(1:round(fitRange*length(msd)));

    tau   = 1:length(tofit);
    f     = fit(tau(:),tofit(:),'a*x');
    
    g = coeffvalues(f);
    D = g(1)/div;

end