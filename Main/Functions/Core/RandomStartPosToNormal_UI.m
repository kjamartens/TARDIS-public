%% Function to obtain a random position from a XYZ-XYZ range
% In the TARDIS UI, a start value can be initialised with a random value,
% e.g. '1-20'. This function extracts a random value out of that range.
%---------------------------------------------------------
% Required inputs
% settingsTARDIS:   settings structure of TARDIS

% Outputs
% settingsTARDIS:   settings structure of TARDIS
% start_1pop:       start values of 1pop fitting
% start_2pop:       start values of 2pop fitting
% start_aDDA:       start values of aDDA fitting
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [settingsTARDIS, start_1pop, start_2pop, start_aDDA] = RandomStartPosToNormal_UI(settingsTARDIS)

%Initialise values if they're not random
start_1pop = settingsTARDIS.start_1pop;
start_2pop = settingsTARDIS.start_2pop;
start_aDDA = settingsTARDIS.start_aDDA;

%Now check in case they are random
if contains(settingsTARDIS.start_1pop,'-')
    %Split on dash - this is between extremes of random initialisation
    t = split(settingsTARDIS.start_1pop(2:end-1),'-');
    settingsTARDIS.start_1pop = ['[' num2str(round(rand()*(str2double(t{2})-str2double(t{1}))+str2double(t{1}),5)) ']'];
    start_1pop = settingsTARDIS.start_1pop;
end

if contains(settingsTARDIS.start_2pop,'-')
    %Create a new running variable with text, where we insert the random
    %pos
    newstart_2pop = '[';
    %Split on space - this is between variables
    s = split(settingsTARDIS.start_2pop(2:end-1),' ');
    for ss = 1:size(s,1)
        if contains(s{ss},'-')
            %Split on dash - this is between extremes of random initialisation
            t = split(s{ss},'-');
            newstart_2pop = [newstart_2pop ' ' num2str(round(rand()*(str2double(t{2})-str2double(t{1}))+str2double(t{1}),5))];
        else
            %If this variable has no dash, just keep it as-is
            newstart_2pop = [newstart_2pop ' ' s{ss}];
        end
    end
    settingsTARDIS.start_2pop = [newstart_2pop ']'];
    start_2pop = settingsTARDIS.start_2pop;
end

if contains(settingsTARDIS.start_aDDA,'-')
    %Create a new running variable with text, where we insert the random
    %pos
    newstart_aDDA = '[';
    %Split on space - this is between variables
    s = split(settingsTARDIS.start_aDDA(2:end-1),' ');
    for ss = 1:size(s,1)
        if contains(s{ss},'-')
            %Split on dash - this is between extremes of random initialisation
            t = split(s{ss},'-');
            newstart_aDDA = [newstart_aDDA ' ' num2str(round(rand()*(str2double(t{2})-str2double(t{1}))+str2double(t{1}),5))];
        else
            %If this variable has no dash, just keep it as-is
            newstart_aDDA = [newstart_aDDA ' ' s{ss}];
        end
    end
    settingsTARDIS.start_aDDA = [newstart_aDDA ']'];
    start_aDDA = settingsTARDIS.start_aDDA;
end