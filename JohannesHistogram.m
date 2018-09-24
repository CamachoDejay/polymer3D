% The purpose of this program is to: 
clear
close all
clc

%% Load block
% Load first file
tmp = load('N:\Project Z-stacks\New stacks\1K_0,25mg_30 nm_sample 2\TIF_1K_0.25 mg_30 nm_sample2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps.mat');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
Table05 = tmp.(names{1});
clear tmp names
% Load second file
tmp = load('N:\Project Z-stacks\New stacks\1K_0,5mg_30 nm_sample 2\tif_1k_0.5 mg_30 nm_Sample 2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps.mat');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
Table075 = tmp.(names{1});
clear tmp names
% Loqd third file
tmp = load('N:\Project Z-stacks\New stacks\1K_1 mg_30 nm_sample 2\tif_1K_1 mg_30 nm_Sample 2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps.mat');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
Table1 = tmp.(names{1});
clear tmp names
%% USER INPUT
nBin = 20;

%% Output block

% plotting all
table05 = [];
table075 = [];
table1 = [];
for i = 1 : length(Table05)
    
    table05  = [table05 ; Table05(i).Data.Area];
    table075 = [table075 ;  Table075(i).Data.Area];
    table1 = [table1 ;  Table1(i).Data.Area];
    
end

% ask for area an put that into totalV
totalV = table05;
% make sure is a col vector
totalV = totalV(:);
% calculate CDF
[midBin,freq,lFreq] = Plotting.lnbin(totalV,nBin);
hist.x = midBin;
hist.y = lFreq/sum(lFreq);
% store CCDF into p05
p05_hist = hist;

% qsk for area an put that into totalV
totalV = table075;
% make sure is a col vector
totalV = totalV(:);
% calculate CDF
[midBin,freq,lFreq] = Plotting.lnbin(totalV,nBin);
hist.x = midBin;
hist.y = lFreq/sum(lFreq);
% store CCDF into p05
p075_hist = hist;

% qsk for area an put that into totalV
totalV = table1;
% make sure is a col vector
totalV = totalV(:);
% calculate CDF
[midBin,freq,lFreq] = Plotting.lnbin(totalV,nBin);
hist.x = midBin;
hist.y = lFreq/sum(lFreq);
% store CCDF into p05
p1_hist = hist;
clear CDF CCDF totalV

%%
figure()

hist = p05_hist;
plot(hist.x,hist.y,'g')
hold on

hist = p075_hist;
plot(hist.x,hist.y,'b')

hist = p1_hist;
plot(hist.x,hist.y,'r')
hold off

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'3 concentrations','hydrogel 1K'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
% ylim([1e-6 1])
% xlim([0.1 1e5])
legend({'0.5 mg/ml','0.75 mg/ml','1 mg/ml'})

%% Condition 1
figure()
hold on
for i = 1 :length(Table05)
    tmp = Table05(i).Data.Area;

    % ask for area an put that into totalV
    totalV = tmp;
    % make sure is a col vector
    totalV = totalV(:);
    % calculate CDF
    [midBin,freq] = Plotting.lnbin(totalV,nBin);
    hist.x = midBin;
    hist.y = freq;
    plot(hist.x,hist.y)
    
end

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','1K PIC_sample 2_0.25 mg/ml'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
% ylim([1e-6 1])
% xlim([0.1 1e5])
legend()

%% Condition 2
figure()
hold on
for i = 1 :length(Table075)
    tmp = Table075(i).Data.Area;

    % ask for area an put that into totalV
    totalV = tmp;
    % make sure is a col vector
    totalV = totalV(:);
    % calculate CDF
    [midBin,freq] = Plotting.lnbin(totalV,nBin);
    hist.x = midBin;
    hist.y = freq;
    
    plot(hist.x,hist.y)
    
end

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','1K PIC_sample 2_0.5 mg/ml'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
% ylim([1e-6 1])
% xlim([0.1 1e5])
legend()

%% Condition 3
figure()
hold on
for i = 1 :length(Table1)
    tmp = Table1(i).Data.Area;

    % ask for area an put that into totalV
    totalV = tmp;
    % make sure is a col vector
    totalV = totalV(:);
    % calculate CDF
     [midBin,freq] = Plotting.lnbin(totalV,nBin);
    hist.x = midBin;
    hist.y = freq;
    
    plot(hist.x,hist.y)
    
end

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','1K PIC_sample 2_1 mg/ml'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
% ylim([1e-6 1])
% xlim([0.1 1e5])
legend()



%%

% figure(2)
% subplot(3,1,1)
% histogram(pores05.Solidity,0.3:0.01:0.99,'normalization','probability')
% title({'Solidity of the pores',' PIC hydrogel (1K; 0.75 mg/ml) '})
% xlabel('Solidity')
% ylabel('normalised probability')
% subplot(3,1,2)
% histogram(pores07.Solidity,0.3:0.01:0.99,'normalization','probability')
% subplot(3,1,3)
% histogram(pores10.Solidity,0.3:0.01:0.99,'normalization','probability')
% 
% 
% figure(3)
% subplot(3,1,1)
% histogram(pores05.Eccentricity,0.1:0.01:0.99,'normalization','probability')
% subplot(3,1,2)
% histogram(pores07.Eccentricity,0.1:0.01:0.99,'normalization','probability')
% subplot(3,1,3)
% histogram(pores10.Eccentricity,0.1:0.01:0.99,'normalization','probability')