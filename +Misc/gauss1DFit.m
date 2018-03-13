function [FitPar,Fit]=gauss1DFit(A,domain)
A=A(:);
domain =domain(:);
[val,ind] = max(A);

%                   Sigma                       mu             A        y0           
lb        = [0.1                             -domain(end)      0        0];
ub        = [5                                domain(end)     3*val      val];
initguess = [0.6                               domain(ind)     val-min(A) min(A)];

FitPar=lsqcurvefit(@Misc.gauss1D,initguess,domain,A,lb,ub);

Fit=FitPar(3)*exp(-((domain-FitPar(2))./(sqrt(2).*FitPar(1))).^2)+FitPar(4);
