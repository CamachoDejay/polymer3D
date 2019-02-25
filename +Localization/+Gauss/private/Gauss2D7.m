function F = Gauss2D7(x,data)
X = data(:,:,1);
Y = data(:,:,2);

F = x(1)*exp(-(((X-x(5)).^2)./(2*x(2).^2)+((Y-x(6)).^2)./(2*x(3)^2)))...
            +x(1)*exp(-(((X-x(7)).^2)./(2*x(2).^2)+((Y-x(8)).^2)./(2*x(3)^2)))...
            +x(1)*exp(-(((X-x(9)).^2)./(2*x(2).^2)+((Y-x(10)).^2)./(2*x(3)^2)))...
            +x(1)*exp(-(((X-x(11)).^2)./(2*x(2).^2)+((Y-x(12)).^2)./(2*x(3)^2)))...
            +x(1)*exp(-(((X-x(13)).^2)./(2*x(2).^2)+((Y-x(14)).^2)./(2*x(3)^2)))...
            +x(1)*exp(-(((X-x(13)).^2)./(2*x(2).^2)+((Y-x(14)).^2)./(2*x(3)^2)))...
            +x(1)*exp(-(((X-x(15)).^2)./(2*x(2).^2)+((Y-x(16)).^2)./(2*x(3)^2)))+x(4);