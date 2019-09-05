%Fitting
function [gpar,resnorm,res] = MultipleGFit3D(A,X0,Y0,Z0,domain,NbG,width)
maxNFit = 5;
if NbG > maxNFit
    error('too many PSF are asked, maximum is 8 for fitting')
end
if length(X0) < maxNFit
    X0(length(X0)+1:maxNFit) = 0;
    Y0(length(Y0)+1:maxNFit) = 0;
    Z0(length(Z0)+1:maxNFit) = 0;
end

if width.xy >0
    
    wguessXY = width.xy;
    lwXY = width.xy;
    uwXY = width.xy;
else
    wguessXY = 3;
    lwXY = 1;
    uwXY = 6;
   
end

if width.z >0
    wguessZ = width.z;
    lwZ = width.z;
    uwZ = width.z;
    
else
    
    wguessZ = 2;
    lwZ = 0.5;
    uwZ = 8;
    
end

bkgub = 0.9  *  max(A(:));
bkglb = 0.95 * min(A(:));
Bkg0 = median(A(:));
% PSFMicroscope=350; %nm
% PSF=PSFMicroscope/160/2.355; %pixel

curvefitoptions = optimset('Display','off');
% LOWER BOUND FOR FITTING PARAMETER
%[peak Int,     sigma X,    sigma Y,        Bkg,            mean X,         mean Y
lb        = [0.5*(max(A(:))-min(A(:)))              lwXY      lwZ         bkglb...
       min(domain(:))        min(domain(:))     min(domain(:))        min(domain(:))...
       min(domain(:))        min(domain(:))     min(domain(:))        min(domain(:))...
       min(domain(:))        min(domain(:))     min(domain(:))        min(domain(:))...
       min(domain(:))        min(domain(:))     min(domain(:))]; 
% UPPER BOUND FOR FITTING PARAMETER
ub        = [1.2*(max(A(:)))  uwXY      uwZ        bkgub ...
    max(domain(:))  max(domain(:))  max(domain(:))    max(domain(:))...
    max(domain(:))  max(domain(:))  max(domain(:))    max(domain(:))...
    max(domain(:))  max(domain(:))  max(domain(:))    max(domain(:))...
    max(domain(:))  max(domain(:))  max(domain(:))];

lb = lb(1,1:4+(NbG)*3);
ub = ub(1,1:4+(NbG)*3);

%INITIAL GUESS
initguess = [max(A(:))-min(A(:))     wguessXY      wguessZ          Bkg0  ...
            X0(1)          Y0(1)    Z0(1) ...
            X0(2)          Y0(2)    Z0(2) ...
            X0(3)          Y0(3)    Z0(3) ... 
            X0(4)          Y0(4)    Z0(4) ...
            X0(5)          Y0(5)    Z0(5) ];
%             X0(5)          Y0(5)    X0(6)      Y0(6) ...
%             X0(7)          Y0(7)    X0(8)      Y0(8)];
        
initguess = initguess(1,1:4+(NbG)*3);

switch NbG
    case 1
           
        [gpar,resnorm,res] = lsqcurvefit(@Gauss3D,...
            initguess,domain,A,lb,ub,curvefitoptions);
        
    case 2
     
        [gpar,resnorm,res] = lsqcurvefit(@Gauss3D2,initguess,...
            domain,A,lb,ub,curvefitoptions);
        
    case 3
    
       [gpar,resnorm,res] = lsqcurvefit(@Gauss3D3,initguess,...
           domain,A,lb,ub,curvefitoptions);
        
    case 4
   
        [gpar,resnorm,res] = lsqcurvefit(@Gauss3D4,initguess,...
            domain,A,lb,ub,curvefitoptions);
        
    case 5
     
        [gpar,resnorm,res] = lsqcurvefit(@Gauss3D5,initguess,...
            domain,A,lb,ub,curvefitoptions);
        
%     case 6
%                 
%         [gpar,resnorm,res] = lsqcurvefit(@Gauss2D6,initguess,...
%             domain,A,lb,ub,curvefitoptions);
%     case 7
%          [gpar,resnorm,res] = lsqcurvefit(@Gauss2D7,initguess,...
%             domain,A,lb,ub,curvefitoptions);
%     case 8
%          [gpar,resnorm,res] = lsqcurvefit(@Gauss2D8,initguess,...
%             domain,A,lb,ub,curvefitoptions);
    otherwise
        error('The current version of the program cannot fit more than 8 2D gaussian at a time ');
end

end