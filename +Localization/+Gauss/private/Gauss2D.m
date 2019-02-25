%Generate a gaussian

function F = Gauss2D(x,data)
X = data(:,:,1);
Y = data(:,:,2);

F = x(1)*exp(-(((X-x(5)).^2)./(2*x(2).^2)+((Y-x(6)).^2)./(2*x(3)^2)))+x(4);
