%% Custom PDF fitting routine - 2 populations with different diffusion and different drifts
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
function output = CustomPDF_Drift_2pop(xdata,parameters)
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
D_first = iterable_parameters(1)/1e12; %Difussion coefficient of the drifting-population
drift_first = iterable_parameters(2)/1e9*tau_value; %This is in nm 
D_second = iterable_parameters(3)/1e12; %Diffusion coefficient of the second populatin
drift_second = iterable_parameters(4)/1e9*tau_value; %This is in nm 
fraction = iterable_parameters(5); %Fraction of all that is drifting (0-1)

loc_unc = settings_parameters.loc_unc;
strobo_frame_time = settings_parameters.strobo_frame_time*tau_value;
frame_time = settings_parameters.frame_time*tau_value;

%Correcting D for loc_unc and strobo frame time (Berglund 2010)
Dcorr_first = D_first+loc_unc^2/frame_time+D_first*(strobo_frame_time/frame_time)*(1/6);
%Obtainting the corrected sigma value
sigma_corr_first = sqrt(2*Dcorr_first*frame_time);
%Obtaining the Rician distribution based on the sigma and drift-based
%offset
Rician_corr_first = besseli(0,xdata*drift_first./sigma_corr_first^2).*(xdata./sigma_corr_first.^2).*...
    exp(-((xdata.^2+drift_first.^2)/(2*sigma_corr_first^2)));

%Correcting D for loc_unc and strobo frame time (Berglund 2010)
Dcorr_second = D_second+loc_unc^2/frame_time+D_second*(strobo_frame_time/frame_time)*(1/6);
%Obtainting the corrected sigma value
sigma_corr_second = sqrt(2*Dcorr_second*frame_time);
%Obtaining the Rician distribution based on the sigma and drift-based
%offset
Rician_corr_second = besseli(0,xdata*drift_second./sigma_corr_second^2).*(xdata./sigma_corr_second.^2).*...
    exp(-((xdata.^2+drift_second.^2)/(2*sigma_corr_second^2)));

output = (Rician_corr_first*fraction+Rician_corr_second*(1-fraction));
end
