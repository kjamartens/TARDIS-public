function TARDIS_analyticalSettings = GenerateTARDISAnalyticalSettings()
%Values are reasonable given a FoV of ~ 500x500 px, 1 px = 100 nm --> AU = 50x50 um.
%% Assessment parameters
TARDIS_analyticalSettings.lambda_bleach = 3; %in frames
TARDIS_analyticalSettings.p_nonLocalized = 0.1; %Chance of a frame blinking on a trajectory (between 0,1)
TARDIS_analyticalSettings.D = 1*1/50; %Diffusion coefficient in AU2/s - 1 um2/s at 50x50 AU is 1/50 AU2/s
TARDIS_analyticalSettings.sigma = 1/2000; %Localization precision in AU - 1/2000 at 50x50 um FoV is 25 nm
TARDIS_analyticalSettings.dt = 0.01; %Frame time in seconds
TARDIS_analyticalSettings.st = 0.01; %Illumination time in seconds

TARDIS_analyticalSettings.d_traj = 2; %Average starting trajectories per frame
TARDIS_analyticalSettings.d_spurious = 5; %Average spurious localizations per frame

TARDIS_analyticalSettings.lambda_bright = 1e99; %Half-time (in frames) fluorophore stays on before it blinks
TARDIS_analyticalSettings.lambda_dark = 1e-99; %Half-time (in frames) fluorophore stays off before it becomes bright again
%% Performance parameters
%Tau-range to calculate/show over, seperately for distribution and
%contribution
TARDIS_analyticalSettings.tauRange_distribution = [1:3];
TARDIS_analyticalSettings.tauRange_contribution = [1:20];
TARDIS_analyticalSettings.movLength = 100; %length of movie in frames
TARDIS_analyticalSettings.dx = 1e-3; %distance steps to take
TARDIS_analyticalSettings.figureCreation = 1; %Whether to create a figure

%% Unchanged parameters
TARDIS_analyticalSettings.dim = 2; %dimensionality
TARDIS_analyticalSettings.xp = [0:TARDIS_analyticalSettings.dx:sqrt(2)]; %x-list to calculate for, in AU
end