function [alpha] = detDiffTypeAlpha(msd,expTime)

    t = (1:length(msd))*expTime;
    
    
    toFit = log(msd);
    
    fitPar = fit(log(t(:)),toFit(:),'a*x+b');
    
    alpha = fitPar.a;
    




end