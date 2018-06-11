clear
close all
clc

load('Z:\Leica TCS SP8X\Leica PC\Johannes\2018\24052018\tif\cf tamra\poreProps.mat')
pores05 = allData;
load('Z:\Leica TCS SP8X\Leica PC\Johannes\2018\24052018\tif\cf atto\poreProps.mat')
pores07 = allData;
load('Z:\Leica TCS SP8X\Leica PC\Johannes\2018\24052018\tif\sted tamra\poreProps.mat')
pores10 = allData;
load('Z:\Leica TCS SP8X\Leica PC\Johannes\2018\24052018\tif\sted atto\poreProps.mat')
pores15 = allData;

clear allData;
%% Plotting
figure(1)
totalV = pores05.Area;
totalV = totalV(:);
[CDF,CCDF] = Plotting.getCDF(totalV);
p05.CCDF = CCDF;

totalV = pores07.Area;
totalV = totalV(:);
[CDF,CCDF] = Plotting.getCDF(totalV);
p07.CCDF = CCDF;

totalV = pores10.Area;
totalV = totalV(:);
[CDF,CCDF] = Plotting.getCDF(totalV);
p10.CCDF = CCDF;

totalV = pores15.Area;
totalV = totalV(:);
[CDF,CCDF] = Plotting.getCDF(totalV);
p15.CCDF = CCDF;

clear CDF CCDF totalV
figure()
CCDF = p05.CCDF;
plot(CCDF.x,CCDF.y,'g')
hold on
CCDF = p07.CCDF;
plot(CCDF.x,CCDF.y,'b')

CCDF = p10.CCDF;
plot(CCDF.x,CCDF.y,'r')

CCDF = p15.CCDF;
plot(CCDF.x,CCDF.y,'y')
hold off

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','2K cf vs sted for 0.5'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
ylim([1e-6 1])
xlim([0.1 1e5])
legend({'CF_tamra','CF_atto 647N','STED_tamra','STED_atto 647N'})
