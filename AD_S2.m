% The purpose of this program is to: 
clear
close all
clc

%% Load block
tmp = load('N:\Project Z-stacks\New stacks\5K_0.25mg_30 nm_Sample 2\tif_5K_0.25mg_30 nm_Sample 2\correct stack\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
Table05 = tmp.(names{1});
clear tmp names
% Load second file
tmp = load('N:\Project Z-stacks\New stacks\5k_0.5 mg 30 nm_Sample 2\tif_5K_0.5mg_30 nm_Sample 2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
Table075 = tmp.(names{1});
clear tmp names
% Loqd third file
tmp = load('N:\Project Z-stacks\New stacks\5K_1 mg_30nm_Sample 2\tif_5K_1 mg_30nm_Sample 2\SegmentedStacks\PoreSize-ResultsAdaptive-poreProps');
names = fieldnames(tmp);
assert(length(names)==1,'Unexpected, WTF')
Table1 = tmp.(names{1});
clear tmp names


%% New Plotting

% plotting all
table05 = [];
table075 = [];
table1 = [];
for i = 1 : length(Table05)
    
    table05  = [table05 ; Table05(i).Data.Area];
    if i<=length(Table075)
        table075 = [table075 ;  Table075(i).Data.Area];
    end
    if i<=length(Table1)
        table1 = [table1 ;  Table1(i).Data.Area];
    end
    
end

% ask for area an put that into totalV
totalV = table05;
% make sure is a col vector
totalV = totalV(:);
% calculate CDF
[CDF,CCDF] = Plotting.getCDF(totalV);
% store CCDF into p05
p05_CCDF = CCDF;

% qsk for area an put that into totalV
totalV = table075;
% make sure is a col vector
totalV = totalV(:);
% calculate CDF
[CDF,CCDF] = Plotting.getCDF(totalV);
% store CCDF into p05
p075_CCDF = CCDF;

% qsk for area an put that into totalV
totalV = table1;
% make sure is a col vector
totalV = totalV(:);
% calculate CDF
[CDF,CCDF] = Plotting.getCDF(totalV);
% store CCDF into p05
p1_CCDF = CCDF;
clear CDF CCDF totalV

%%
figure()

CCDF = p05_CCDF;
plot(CCDF.x,CCDF.y,'g')
hold on

CCDF = p075_CCDF;
plot(CCDF.x,CCDF.y,'b')

CCDF = p1_CCDF;
plot(CCDF.x,CCDF.y,'r')
hold off

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','5K PIC hydrogel for different concentrations (adaptive) S2'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
ylim([1e-6 1])
xlim([0.1 1e5])
legend({'0.25 mg/ml','0.5 mg/ml','1 mg/ml'})

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
    [CDF,CCDF] = Plotting.getCDF(totalV);
    
    plot(CCDF.x,CCDF.y)
    
end

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','5K PIC hydrogel 0.25 mg/ml adaptive sample 2'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
ylim([1e-6 1])
xlim([0.1 1e5])
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
    [CDF,CCDF] = Plotting.getCDF(totalV);
    
    plot(CCDF.x,CCDF.y)
    
end

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','5K PIC hydrogel 0.5 mg/ml adaptive sample 2'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
ylim([1e-6 1])
xlim([0.1 1e5])
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
    [CDF,CCDF] = Plotting.getCDF(totalV);
    
    plot(CCDF.x,CCDF.y)
    
end

a = gca;
a.XScale = 'log';
a.YScale = 'log';
title({'CCDF-plot for area of pores','5K PIC hydrogel 1 mg/ml adaptive sample 2'})
xlabel('Pore size (micrometer^2)')
ylabel('CCDF - prob [0-1]')
a.FontSize = 14;
ylim([1e-6 1])
xlim([0.1 1e5])
legend()