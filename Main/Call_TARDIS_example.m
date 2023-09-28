%Example for using TARDIS from command line
clear all; 
close all;
%% Load data
%Loaded data is a 'pos' array with frame, x, y position (x, y, in meter)
load('.\Example_data\Single_Population_1um2s_lowDensity.mat','pos');

%% Initialise settings
settingsTARDIS = GenerateTARDISsettings();
%User-define the most important settings
settingsTARDIS.maxdist = 5e-6; %Maximum dist (m)
settingsTARDIS.startpointBG = 3.2e-6;
settingsTARDIS.populations = 1; %number of populations to fit
settingsTARDIS.loc_unc = 30*1e-9; %Localization uncertainty in m
settingsTARDIS.fitWithBleach = 1; %Bleach fit or not
settingsTARDIS.performestimationfit = 0; %Estimation fit or not
settingsTARDIS.linorlogvis = 'lin'; %Linear or logartihmic fit+visualisation

%% Determine the BG frame via the Wilcoxon test and change the settings accordingly
%This is done for 36 frame-delays, repeating 5 times, no visualisation
settingsTARDIS.frame_dist_BG_start = BGLengthDetermination_Wilcoxon(pos,36,5,settingsTARDIS.maxdist,0);
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
disp(parameters);

%% Visualise the fit - total fit visualisation is done via the FullFit (FF) visualisation:
% In case you want to repeat the visualisation, you'd need the
% Visualisation_FF_outputCell structure like this:
Visualisation_FF([],[],Visualisation_FF_outputCell);

%% Visualise the JD histogram visualisation
%Create JD histogram visualisation
%Note that this basically just histograms the info in JDonlydata
extrainfo.Data.settingsURDA = settingsTARDIS;
Visualisation_JDs('','',JDonlydata{1},extrainfo);

%% Run TARDIS GUI
%Alternatively, you can open the TARDIS GUI with this:
TARDIS_app;

%% You could also re-create output images with the stored MAT file like this:
%The TARDIS APP stores .mat files as output in a TARDIS_Results Folder.
%Load this data first!
%First we load the .mat TARDIS GUI output
%And we plot image like this:
Visualisation_FF([],[],appinfo.Visualisation_FF_outputCell);