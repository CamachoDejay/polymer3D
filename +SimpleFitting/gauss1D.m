function [FitPar,Fit]=gauss1D(A,domain,guess)

switch nargin
    case 2 
        sigGuess = (abs(domain(2))-abs(domain(1)))*3;
        [val,idx] = max(A);
        muGuess = domain(idx);
    case 3
        sigGuess = guess.sig;
        muGuess  = guess.mu;
    otherwise
        error('wrong number of arguments');
end

A=A(:);
domain =domain(:);
[val,ind] = max(A);

%                   Sigma                       mu             A        y0           
lb        = [abs(domain(1)-domain(2))/10     min(domain)-0.1*abs(max(domain))     0        0];
ub        = [abs(domain(1)-domain(2))*20     max(domain)+0.1*abs(max(domain))     3*val      val];
initguess = [      sigGuess                  muGuess                          val-min(A) min(A)];

FitPar=lsqcurvefit(@SimpleFitting.gaussian,initguess,domain,A,lb,ub);

Fit=FitPar(3)*exp(-((domain-FitPar(2))./(sqrt(2).*FitPar(1))).^2)+FitPar(4);
