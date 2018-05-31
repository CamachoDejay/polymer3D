%testlocalization: The aim of the code is to test all the function linked to
%localization and fitting that are used in the polymer3D project.

%   The user will be asked which part of the code he want to test (several
%   option and all will be proposed). Then the code will run the chosen
%   test function(s) that are located in +test\ folder.
clear;
close all;
clc;
%% Radio button for option choice
fig = uifigure('Name','Please choose the unit to be tested','Position',[500 250 600 300]);
bg = uibuttongroup(fig,'Position',[50 50 500 200]);   

rb1 = uiradiobutton(bg,'Position',[10 180 500 20]);
rb1.Text = 'Test Rough Localization';
rb2 = uiradiobutton(bg,'Position',[10 100 500 20]);
rb2.Text = 'Test subdiffraction Localization';
rb3 = uiradiobutton(bg,'Position',[10 20 500 20]);
rb3.Text = 'Test full procedure (from rough localization on simulated images to fitting)';

option = [rb1,rb2,rb3];
b4 = uibutton(bg,'push', 'Position',[440,20,40,20],'ButtonPushedFcn', @(btn,event) executeTest(btn, fig, option));
b4.Text =  'Go!';
%% Run the chosen test
function executeTest(~,fig, option)
test2Run = [option(1).Value,option(2).Value,option(3).Value];
[~,ind] = find(test2Run);
close(fig);

switch ind
    case 1
        disp('Starting the test of GLRT localization---------------------');
        Test.GLRTLocalization;
        disp('--------------------- Test Done ! ---------------------');
    case 2
        disp('Starting the test of subdiffraction localization-----------');
        Test.SRLocalization;
        disp('--------------------- Test Done ! ---------------------');
    case 3
        disp('Starting the full localization test------------------------');
        disp('--------------------- Test Done ! ---------------------');
end

end