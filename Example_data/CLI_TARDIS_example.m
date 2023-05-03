%Example for using TARDIS from command line
clear all; 
close all;
%% Load data
%Loaded data is a 'pos' array with frame, x, y position (x, y, in meter)
load(['Single_Population_1um2s_lowDensity.mat']);

%% Initialise settings
%For all settings, read GenerateTARDISsettings.
settingsTARDIS = GenerateTARDISsettings();
%User-define the most important settings
settingsTARDIS.maxdist = 3e-6; %Maximum dist (m)
settingsTARDIS.startpointBG = 2.4e-6; %starting position of background-only distribution (m)
settingsTARDIS.populations = 1; %number of populations to fit
settingsTARDIS.loc_unc = 30*1e-9; %Localization uncertainty in m
settingsTARDIS.fitWithBleach = 1; %Bleach fit or not
settingsTARDIS.performestimationfit = 0; %Estimation fit or not
settingsTARDIS.linorlogvis = 'log'; %Linear or logartihmic fit+visualisation

%% Determine the BG frame via the Wilcoxon test and change the settings accordingly
% settingsTARDIS.frame_dist_BG_start = BGLengthDetermination_Wilcoxon(pos,36,5,settingsTARDIS.maxdist,0);
settingsTARDIS.frame_dist_BG_start = 30; %Start frame-delay for background-only population. Alternatively, use Wilcoxon as above
% Using 50 different dt's for BG curve determination
settingsTARDIS.frame_dist_BG = [settingsTARDIS.frame_dist_BG_start+1:settingsTARDIS.frame_dist_BG_start+51];

%% Run TARDIS
%Easiest way to run TARDIS is simply like this:
TARDIS(pos,settingsTARDIS);

%% Run TARDIS with parameter output
%TARDIS outputs are ordered like this
[parameters, parametersCI, paramEsts, HOfitCI, tottime, time, bgarr, truthoffsetpartial, ...
    Visualisation_HO_outputCell,Visualisation_FF_outputCell,anaDDAvisInfoHO,anaDDAvisInfoFF,JDonlydata,SWIFTParams] = ...
    TARDIS(pos,settingsTARDIS);

%The fit information in stored in parameters and parametersCI
disp('Found parameters: (Diffusion coefficient, % of links belonging to background, bleach time (s)')
disp(parameters);

%% Visualise the fit - total fit visualisation is done via the FullFit (FF) visualisation:
% In case you want to repeat the visualisation, you'd need the
% Visualisation_FF_outputCell structure like this:
Visualisation_FF([],[],Visualisation_FF_outputCell);

%% Run TARDIS GUI
%Alternatively, you can open the TARDIS GUI via TARDIS_app