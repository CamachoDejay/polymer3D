function [FitPar,Fit]=gauss1D(A,domain)
A=A(:);
domain =domain(:);
[val,ind] = max(A);

%                   Sigma                       mu             A        y0           
lb        = [abs(domain(1)-domain(2))/10       -2*abs(median(domain))     0        0];
ub        = [abs(domain(1)-domain(2))*20     2*abs(median((domain)))     3*val      val];
initguess = [(domain(1)-domain(2))          domain(ind)     val-min(A) min(A)];

FitPar=lsqcurvefit(@SimpleFitting.gaussian,initguess,domain,A,lb,ub);

Fit=FitPar(3)*exp(-((domain-FitPar(2))./(sqrt(2).*FitPar(1))).^2)+FitPar(4);
