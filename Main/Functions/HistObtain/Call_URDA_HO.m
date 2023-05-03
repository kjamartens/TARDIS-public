load('C:\Users\Koen Martens\Documents\WUR info - data\Backup H Drive\MATLAB_BIPNAS\RelativeDistance_trackingAlternative\SimulatedTracks\OnePop_anaDDA\dens_2_kon_0_koff_0_Dfree_100.mat');

% clc
clear settingsURDA
settingsURDA.frame_dist_BG = [9:1:30];
settingsURDA.startpointBG = 1.5e-6;
settingsURDA.maxdist = 3e-6;
settingsURDA.dt_arr = [1 2 3];
settingsURDA.vis = 1;
settingsURDA.verbose = 1;
settingsURDA.fitvis = 1;
settingsURDA.bgbinningnr = 100;
settingsURDA.linorlogBGsubtract = 'lin';
settingsURDA.linorlogvis = settingsURDA.linorlogBGsubtract;
settingsURDA.minlogpoint = 10^-9;
settingsURDA.frame_time = 0.01;
settingsURDA.loc_unc = 15e-9;
settingsURDA.norm_bins = 10000;
settingsURDA.start = [1e-12]; %Startpoint of fitting
settingsURDA.lb = [1e-16]; %LB of fit
settingsURDA.ub = [30e-12]; %UB of fit
settingsURDA.debug = false;

URDA_HO_MLE_function(poslist,settingsURDA);
