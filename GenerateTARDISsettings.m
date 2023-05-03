function TARDISsettings = GenerateTARDISsettings()
% 'Crucial' settings
TARDISsettings.frame_time = 0.01; %- %Frame time in seconds
TARDISsettings.loc_unc = 15e-9; %- %Localization uncertainty in meters
TARDISsettings.populations = 1; %- %number of populations to be fit. Make this 0 to fit anaDDA!
TARDISsettings.AutoChoosePop = 0; %Automatically choose between 1, 2, or 3 (not yet implemented 3) populations
TARDISsettings.frame_dist_BG = [20:40]; %- % Which dt's are used to get BG (f.e. [9:1:20])
TARDISsettings.dt_arr = [1 2 3]; %- %which dt's are used to get signal
TARDISsettings.performsecondfit = 1; %Perform second (better) MLE fit

% 'Medium' settings
TARDISsettings.performestimationfit = 1; %- %Perform estimation fit
TARDISsettings.maxdist = 3e-6; %- %What the max JD is
TARDISsettings.fitWithBleach = 0; %fit BG/populations with bleach kinetics rather than with arbitrary values
TARDISsettings.fixRatios2pop = 1; %Fit 2pop with fixed ratios (true) or not (false)
TARDISsettings.startpointBG = 1e-6; %- %what JD the BG-only area is starting
TARDISsettings.start_1pop = '[0.1-20]'; %- %Startpoint of fitting between 1e-11 and 1e-13
TARDISsettings.lb_1pop = '[1e-5]'; %LB of fit
TARDISsettings.ub_1pop = '[50]'; %UB of fit
TARDISsettings.start_2pop = '[0.1-20 0.1-20 0.1-10]';%-  %Startpoint of fitting between 1e-11 and 1e-13
TARDISsettings.lb_2pop = '[1e-5 1e-5 1e-5]'; %LB of fit
TARDISsettings.ub_2pop = '[50 50 100]'; %UB of fit
TARDISsettings.start_3pop = '[0.05 1 29 0.9]'; %Startpoint of fitting between 1e-11 and 1e-13
TARDISsettings.lb_3pop = '[1e-5 1e-5 1e-5 0]'; %LB of fit
TARDISsettings.ub_3pop = '[1e4 1e-4 1e4 100]'; %UB of fit
TARDISsettings.start_aDDA = '[10-80 10-80 0.1-10]'; %- %Startpoint of anaDDAfit
TARDISsettings.lb_adda = '[1e-5 1e-5 1e-5]';
TARDISsettings.ub_adda = '[1e4 1e4 1e2]';
TARDISsettings.stroboFrameTime = 0; %Stroboscopic illumination in seconds. Set to 0 for no strobo-illumination (i.e. ignored), -1 to full-frame-illumination, or a specific time.
TARDISsettings.freefit_locunc = 0; %EXPERIMENTAL free fitting of localization uncertainty

% 'Low importance' settings
TARDISsettings.bgbinningnr = 100; %- %nr of BG bins
TARDISsettings.linorlogBGsubtract = 'log'; %'lin' or 'log'
TARDISsettings.linorlogvis = TARDISsettings.linorlogBGsubtract; %- ;
TARDISsettings.minlogpoint = 10^(-8.5); %starting point if using log BG subtract - used to be -8.5
TARDISsettings.norm_bins = 50000; %used to be 10000 %Nr of bins used for normalization - everything > 1000 seems fine
TARDISsettings.createJDsonly = 0; %Only create a JD list rather than try to fit something
TARDISsettings.noiseDensity = 0; %density (loc/um2) of noisy (i.e. non-track) particles. Used for determining blink chance

% 'Verbose'-related settings
TARDISsettings.verbose = 1; 
TARDISsettings.verboseReal = TARDISsettings.verbose; %actual verbose of fitting etc
TARDISsettings.debug = 0; %Debug info
TARDISsettings.vis = TARDISsettings.debug; %Visibility of BG subtracting results
TARDISsettings.fitvisHO = TARDISsettings.debug; %Visibility of fit results
TARDISsettings.visualisationMLEIntFit = 1; %Visualisation of MLE fit at the end
TARDISsettings.StoreSWIFTparameters = 0; %Store a JSON for SWIFT settings
TARDISsettings.verboseMLEIntFit = 0;
TARDISsettings.callfromUI = 0;

%% aDDA settings
TARDISsettings.inputAnaDDA = Generateinputfile;
TARDISsettings.inputAnaDDA.numberofspecies = 1;
TARDISsettings.inputAnaDDA.confinement = 1; %in app
TARDISsettings.inputAnaDDA.radiusofcell = 0.5; %in app
TARDISsettings.inputAnaDDA.lengthcell = 3; %in app

TARDISsettings.inputAnaDDA.upperDfree = 10; % in um2/s
TARDISsettings.inputAnaDDA.compensatetracking = 0;
TARDISsettings.inputAnaDDA.fixedparameters(:) = -1;
TARDISsettings.inputAnaDDA.framerange = 1; %hardcoded and should always be this way
TARDISsettings.inputAnaDDA.precision = 2^16; %Normally 50000

% TARDISsettings.useAIC = 0; %Using AIC to elucidate the number of populations
% TARDISsettings.callfromUI = 0; %Set to the UI app (callfromUI = app) if called from UI - usefull for text output;
end