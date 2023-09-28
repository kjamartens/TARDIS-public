%Example on how to use TARDIS in a script
TARDISsettings = GenerateTARDISsettings();

%Take a second look at the most crucial settings
TARDISsettings.frame_time = 0.01; %- %Frame time in seconds
TARDISsettings.loc_unc = 15e-9; %- %Localization uncertainty in meters
TARDISsettings.populations = 1; %- %number of populations to be fit. Make this 0 to fit anaDDA!
TARDISsettings.AutoChoosePop = 0; %Automatically choose between 1, 2, or 3 (not yet implemented 3) populations
TARDISsettings.frame_dist_BG = [20:40]; %- % Which dt's are used to get BG (f.e. [9:1:20])
TARDISsettings.dt_arr = [1 2 3]; %- %which dt's are used to get signal
TARDISsettings.performsecondfit = 1; %Perform second (better) MLE fit
TARDISsettings.performestimationfit = 0; %Perform estimation
TARDISsettings.maxdist = 3e-6;
TARDISsettings.visualisationMLEIntFit = 1; %Visualisation of MLE fit at the end

TARDISsettings.debug = 0; %Debug info
TARDISsettings.vis = TARDISsettings.debug; %Visibility of BG subtracting results
TARDISsettings.fitvisHO = TARDISsettings.debug; %Visibility of fit results
TARDISsettings.bgbinningnr = 100;

%Specify file location
fileloc = '.\Example_data\Single_Population_1um2s_lowDensity.mat';
%Load the data and store as poslist
t = load(fileloc,'pos');
poslist = t.pos;

%Run TARDIS
[parameters, parametersCI, paramEsts, HOfitCI, tottime, time, bgarr, truthoffsetpartial, Visualisation_HO_outputCell,Visualisation_FF_outputCell,anaDDAvisInfoHO,anaDDAvisInfoFF,JDonlydata] = ...
    TARDIS(poslist,TARDISsettings);

%% Data visualisation/showcasing
%Create final fit visualisation
Visualisation_FF('','',Visualisation_FF_outputCell);

%Create estimation fit visualisation
extrainfo.Data.settingsURDA = TARDISsettings;
Visualisation_HO_BGsubtract('','',Visualisation_HO_outputCell.bgbinningnr,Visualisation_HO_outputCell.dt_arr,...
    Visualisation_HO_outputCell.JDarrBG,Visualisation_HO_outputCell.maxdist,Visualisation_HO_outputCell.minlogpoint,...
    Visualisation_HO_outputCell.linorlogBGsubtract,Visualisation_HO_outputCell.bgarr,Visualisation_HO_outputCell.JDarrSignalCell,...
    Visualisation_HO_outputCell.signalCurve_interp,Visualisation_HO_outputCell.stepsizearraytrue,Visualisation_HO_outputCell.bgratios,Visualisation_HO_outputCell.startpointBG,...
    extrainfo);

%Create JD histogram visualisation
Visualisation_JDs('','',JDonlydata{1},extrainfo);

            