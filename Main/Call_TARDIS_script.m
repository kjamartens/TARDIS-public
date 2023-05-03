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
fileloc = '\\IFMB-NAS\AG Endesfelder\Data\Koen\Articles\2020_relativeDisplacementTracking\Data_FINAL_202201\1pop_Density\TestCopy\nocutoff_singlePop_densMult10_1_locum2_DfreeMult10_10_um2s.mat';
% fileloc = '\\ifmb-nas\ag endesfelder\Data\Koen\Articles\2020_relativeDisplacementTracking\Data_FINAL_202201\aDDA_Noise_Offset_Bias\aDDANoise_20_20_2_oLoc_0_sLoc_0_it_1_2kframes.mat'
t = load(fileloc,'pos');
poslist = t.pos;
%
[parameters, parametersCI, paramEsts, HOfitCI, tottime, time, bgarr, truthoffsetpartial, Visualisation_HO_outputCell,Visualisation_FF_outputCell,anaDDAvisInfoHO,anaDDAvisInfoFF,JDonlydata] = ...
    URDA_HO_FF_function(poslist,TARDISsettings);

%%
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

            