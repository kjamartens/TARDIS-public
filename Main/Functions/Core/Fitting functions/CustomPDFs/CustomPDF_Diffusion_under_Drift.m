%% Custom PDF fitting routine - Diffusion under drift
% The actual fit called by e.g. HistOnly_MLEFit_estimation or
% pdfSinglePopFunction
%---------------------------------------------------------
% Inputs
% xdata:                All input values of the JD data
% parameters:           A set of parameters - includes the iterable
%                       parameters, settings parameters, and the current
%                       temporal shift
%
% Output
% output                Likelihood of every xdata to belong to this
%                       PDF/params
%---------------------------------------------------------
% Koen J.A. Martens, 2023
%---------------------------------------------------------
function output = CustomPDF_Diffusion_under_Drift(xdata,parameters)
%Get the iterable-parameters out - this is what the PDF should be fitting!
%Note that TARDIS is externally also fitting intra-emitter fraction and
%bleach-time.
iterable_parameters = parameters{1};
%Get settings information
settings_parameters = parameters{2};
%Get the current tau value (temporal shift)
tau_value = parameters{3};

%% Custom PDF area starts here
%Get readable parameters
D = iterable_parameters(1)/1e12;
drift = iterable_parameters(2)/1e9*tau_value; %This is in nm
loc_unc = settings_parameters.loc_unc;
strobo_frame_time = settings_parameters.strobo_frame_time*tau_value;
frame_time = settings_parameters.frame_time*tau_value;

%Correcting D for loc_unc and strobo frame time (Berglund 2010)
Dcorr = D+loc_unc^2/frame_time+D*(strobo_frame_time/frame_time)*(1/6);
%Obtainting the corrected sigma value
sigma_corr = sqrt(2*Dcorr*frame_time);
%Obtaining the Rician distribution based on the sigma and drift-based
%offset
Rician_corr = besseli(0,xdata*drift./sigma_corr^2).*(xdata./sigma_corr.^2).*...
    exp(-((xdata.^2+drift.^2)/(2*sigma_corr^2)));

output = Rician_corr;
end
