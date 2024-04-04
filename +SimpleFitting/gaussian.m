function F=gaussian(x,data)

F = x(3) .* exp(-((data-x(2))./(sqrt(2).*x(1))).^2)+x(4);




% for testing 
% test = randn(10000,1);
% 
% [N,edges] = histcounts(test,100);
% figure
% plot(edges(1:end-1),N);