%% Generate typical TARDIS settings
% Good for start for commandline-based TARDIS fitting
%---------------------------------------------------------
% Required inputs
% none
%
% Output
% settingsTARDIS        TARDIS settings structure
% Last updated: 2021-12-13
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function settingsTARDIS = generateTARDISsettings()
settingsTARDIS.createJDsonly = 0;
settingsTARDIS.maxdist = 3.00e-06;
settingsTARDIS.dt_arr = [1,2,3];
settingsTARDIS.frame_time = 0.01;
settingsTARDIS.frame_dist_BG_start = 12;
settingsTARDIS.frame_dist_BG = [settingsTARDIS.frame_dist_BG_start+1:settingsTARDIS.frame_dist_BG_start+51];
settingsTARDIS.performestimationfit = 0;
settingsTARDIS.verboseReal = 1;
settingsTARDIS.bgbinningnr = 300;
settingsTARDIS.stroboFrameTime = 0;
settingsTARDIS.startpointBG = 2.55e-06;
settingsTARDIS.fitWithBleach = 1;
settingsTARDIS.linorlogvis = 'log';
settingsTARDIS.loc_unc = 4.00e-08;
settingsTARDIS.populations = 1;
settingsTARDIS.start_1pop = '[0.1-10]';
settingsTARDIS.start_2pop = '[.01-100.01-100.2-5]';
settingsTARDIS.start_aDDA = '[2000.12]';
settingsTARDIS.savemaxdata = 0;
settingsTARDIS.linorlogBGsubtract = 'log';
settingsTARDIS.inputAnaDDA = Generateinputfile_KMmodified;
settingsTARDIS.inputAnaDDA.sigmaerror = settingsTARDIS.loc_unc*1e6;
settingsTARDIS.fitvis = 1;
settingsTARDIS.verbose = 1;
settingsTARDIS.minlogpoint = 1e-08;
settingsTARDIS.norm_bins = 10000;
settingsTARDIS.lb_1pop = 1e-05;
settingsTARDIS.ub_1pop = 10000;
settingsTARDIS.lb_2pop = [1e-05,1e-05,0];
settingsTARDIS.ub_2pop = [100,100,100];
settingsTARDIS.start_3pop = [0.05,1,29,0.9];
settingsTARDIS.lb_3pop = [1e-05,1e-05,1e-05,0];
settingsTARDIS.ub_3pop = [10000,10000,10000,100];
settingsTARDIS.start_aDDA = [30,30,2];
settingsTARDIS.lb_aDDA = [1e-5, 1e-5, 1e-5];
settingsTARDIS.ub_aDDA = [1e4, 1e4, 1e2];
settingsTARDIS.freefit_locunc = False;
settingsTARDIS.debug = 0;
settingsTARDIS.performsecondfit = 1;
settingsTARDIS.verboseMLEIntFit = 0;
settingsTARDIS.vis = 0;
settingsTARDIS.fitvisHO = 0;
settingsTARDIS.visualisationMLEIntFit = 0;
settingsTARDIS.useAIC = 0;
settingsTARDIS.callfromUI = 0;
