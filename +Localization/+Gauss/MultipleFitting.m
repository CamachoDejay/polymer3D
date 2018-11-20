%Fitting
function [gpar,resnorm,res] = MultipleFitting(A,X0,Y0,domain,NbG)
maxNFit = 6;
if NbG > maxNFit
    error('too many PSF are asked, maximum is 6 for fitting')
end
if length(X0) < maxNFit
    X0(length(X0)+1:maxNFit) = 0;
    Y0(length(Y0)+1:maxNFit) = 0;
end
radius=(size(A,2)-1)/2;

wguess = 3;
lw = 1;
uw = 4;

bkgub = 2*min(A(:));
bkglb = 0.95*min(A(:));
Bkg0 = median(A(:));
% PSFMicroscope=350; %nm
% PSF=PSFMicroscope/160/2.355; %pixel

curvefitoptions = optimset('Display','off');

%[peak Int,     sigma X,    sigma Y,        Bkg,            mean X,         mean Y
lb        = [0.5*(max(A(:))-min(A(:)))              lw          lw               bkglb...
       min(domain(:))        min(domain(:))     min(domain(:))        min(domain(:))...
       min(domain(:))        min(domain(:))     min(domain(:))        min(domain(:))...
       min(domain(:))        min(domain(:))     min(domain(:))        min(domain(:)) ];

ub        = [1.2*(max(A(:))-min(A(:)))  uw      uw        bkgub ...
    max(domain(:))  max(domain(:))  max(domain(:))    max(domain(:))...
    max(domain(:))  max(domain(:))  max(domain(:))    max(domain(:))...
    max(domain(:))  max(domain(:))  max(domain(:))    max(domain(:))...
    max(domain(:))  max(domain(:))  max(domain(:))    max(domain(:))];

lb = lb(1,1:6+(NbG-1)*2);
ub = ub(1,1:6+(NbG-1)*2);

initguess = [max(A(:))-min(A(:))     wguess      wguess          Bkg0  ...
            X0(1)          Y0(1)    X0(2)      Y0(2) ...
            X0(3)          Y0(3)    X0(4)      Y0(4) ...
            X0(5)          Y0(5)    X0(6)      Y0(6) ];
        
initguess = initguess(1,1:6+(NbG-1)*2);

switch NbG
    case 1
           
        [gpar,resnorm,res] = lsqcurvefit(@Localization.Gauss.Gauss2D,...
            initguess,domain,A,lb,ub,curvefitoptions);
        
    case 2
     
        [gpar,resnorm,res] = lsqcurvefit(@Localization.Gauss.Gauss2D2,initguess,...
            domain,A,lb,ub,curvefitoptions);
        
    case 3
    
       [gpar,resnorm,res] = lsqcurvefit(@Localization.Gauss.Gauss2D3,initguess,...
           domain,A,lb,ub,curvefitoptions);
        
    case 4
   
        [gpar,resnorm,res] = lsqcurvefit(@Localization.Gauss.Gauss2D4,initguess,...
            domain,A,lb,ub,curvefitoptions);
        
    case 5
     
        [gpar,resnorm,res] = lsqcurvefit(@Localization.Gauss.Gauss2D5,initguess,...
            domain,A,lb,ub,curvefitoptions);
        
    case 6
                
        [gpar,resnorm,res] = lsqcurvefit(@Localization.Gauss.Gauss2D6,initguess,...
            domain,A,lb,ub,curvefitoptions);
        
    otherwise
        error('The current version of the program cannot fit more than 6 2D gaussian at a time ');
end

end