%Obtain settings
TARDIS_analyticalSettings = GenerateTARDISAnalyticalSettings();
%Ensure visualisation is on
TARDIS_analyticalSettings.figureCreation = 1;
TARDIS_analyticalSettings.lambda_bright = 3; %Rate (in frames) fluorophore stays on before it blinks
TARDIS_analyticalSettings.lambda_dark = 2; %Blinking rate (in frames)
%Run TARDIS_analytical
[distribution_TARDIS,distribution_interParticle,distribution_intraParticle,...
    fraction_interParticle,xdata,n_linkages] = TARDIS_analytical(TARDIS_analyticalSettings);

%% Alternatively, run TARDIS_analytical_syms, which uses the symbolics
%toolbox
%The TARDIS_syms output is a formula for the TARDIS plot
[TARDIS_syms] = TARDIS_analytical_syms(TARDIS_analyticalSettings);
