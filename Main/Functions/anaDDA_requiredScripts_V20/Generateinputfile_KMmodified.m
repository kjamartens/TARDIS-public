function [input] = Generateinputfile_KMmodified
%% This script generates the input parameters that are used in both fitting analysis and simulations.
% Can be run before running anaDDA to avoid prompt or to modify more advanced settings not modifiable in the prompt. To load this into the
% anaDDA pipeline, first type 'input = Generateinputfile' and subsequently
% 'anaDDA(input)'. 

%% Basic input parameters (Can also be altered via the prompt if anaDDA is run without inputfile)
input.numberofspecies = 1;                                      % Number of (non interacting) species you want to fit (max 3)
input.upperDfree = 10;                                          % Upper estimate of Dfree (um2/s)
input.sigmaerror = 0.03;                                        % Localization error used in both fitting and simulations(um)
input.confinement = false;                                      % Whether confinement boundaries are taken in account in fitting data
input.radiusofcell = 0.5;                                       % The radius of the spherical/rod-shaped cells/boundaries used for both fitting and simulations (um)
input.lengthcell = 3;                                           % The length of the spherical/rod-shaped cells/boundaries used for both fitting and simulations (um)
input.compensatetracking = true;                                % Whether tracking windows are taken in account in fitting data
input.trackingwindow = 300;                                     % length of tracking window (um)
input.fixedparameters = [1 -1 -1 -1; -1 -1 -1 -1;-1 -1 -1 -1];  % Fixed parameters for fitting. Each row is a species, the columns are fraction, koff, kon and Dfree. -1 if not fixed(Default: [1 -1 -1 -1; -1 -1 -1 -1;-1 -1 -1 -1];)

%% Advanced fitting parameters
input.bootstrapping = true;                                     % If you want to include bootstrapping for std dev estimations
input.numberofbootstraps = 8;                                   % Number of bootstraps
input.cyclenumber = 4;                                          % Minimum amount of rounds of MLE before comparing parameter outputs and potential termination
input.lowerstartkoff = 0.05;                                    % Lower bound of koff estimate in MLE optimization (s-1)
input.upperstartkoff = 2000;                                    % Upper bound of koff estimate in MLE optimization (s-1)
input.lowerstartkon = 0.1;                                      % Lower bound of kon estimate in MLE optimization (represented as fraction of koff) 
input.upperstartkon = 10;                                       % Upper bound of kon estimate in MLE optimization (represented as fraction of koff) 
input.precision = 50000;                                        % Number of points in distribution for which anaDDA calculates likelihood directly (instead of interpolation)
input.nofit = false;                                            % Whether fitting algorithm is used (false) or directly use input parameters from simulation to calculate distribution (true) 
input.plotlog = true;                                           % Whether plots are being generated (true)
input.KSstats = true;                                           % Whether Kolmogorov Smirnov test statistic is calculated (true)

%% Advanced simulation parameters
% These parameters only matter in case you use simulation
input.koff1_A = 0;                                             % first off rate for species one (s-1)
input.koff2_A = 10e9;                                           % second off rate for species one (if no second off rate 10e9/much higher than first on rate)
input.kon1_A = 0;                                              % first on rate for species one (s-1)
input.kon2_A = 0.0001;                                          % second off rate for species one (if no second on rate 0.0001/much lower than first on rate)
input.Dfree_A = 1;                                              % free diffusion coefficient for species one (um2/s)

input.fractionB = 0;                                            % first off rate for species one (s-1)
input.koff1_B = 0;                                             % first off rate for species two (s-1)
input.koff2_B = 10e9;                                           % second off rate for species two (s-1)
input.kon1_B = 0;                                              % first on rate for species two (s-1)
input.kon2_B = 0.0001;                                          % second on rate for species two (s-1)
input.Dfree_B= 4.5;                                             % free diffusion coefficient for species two (um2/s)
input.Nparticles = 50;                                       % Number of particles (tracks) run in the simulation
input.distributionNparticles = [0.28723 0.20581 0.14747 0.10567 0.07572 0.05425 0.03887 0.08496]; % Distribution of track lengths for each length ranging from 1 to 8 steps
input.framerange = [1 2 3 4 5 6 7 8];                           % Which track lengths are simulated (number of steps)
input.pixelsize = 1;
input.Steptime = 0.000001; %0.000001 is 1 us                                       % step time of simulation (s) 
input.frametimerange = 0.01;                                    % Range of frame times for which particles are simulated (s)

% %IDEA FOR LATER: Can be easily expanded for the simulation!
% ii = 1;
% input.Species(ii).fracttion = 1; 
% input.Species(ii).koff1 = 30;                                             % first off rate for species one
% input.Species(ii).koff2 = 10e9;                                           % second off rate for species one (if no second off rate 10e9/much higher than first on rate)
% input.Species(ii).kon1 = 30;                                              % first on rate for species one
% input.Species(ii).kon2 = 0.0001;                                          % second off rate for species one (if no second on rate 0.0001/much lower than first on rate)
% input.Species(ii).Dfree = 1;                                              % first off rate for species one
% 
% ii=ii+1
% input.Species(ii).fracttion = 1; 
% input.Species(ii).koff1 = 30;                                             % first off rate for species one
% input.Species(ii).koff2 = 10e9;                                           % second off rate for species one (if no second off rate 10e9/much higher than first on rate)
% input.Species(ii).kon1 = 30;                                              % first on rate for species one
% input.Species(ii).kon2 = 0.0001;                                          % second off rate for species one (if no second on rate 0.0001/much lower than first on rate)
% input.Species(ii).Dfree = 4.5;                                              % first off rate for species one
