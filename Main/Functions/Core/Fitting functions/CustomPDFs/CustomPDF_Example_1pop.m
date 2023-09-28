%% Custom PDF fitting routine - example of 1 population
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
function output = CustomPDF_Example_1pop(xdata,parameters)
%Get the iterable-parameters out - this is what the PDF should be fitting!
%Note that TARDIS is externally also fitting intra-emitter fraction and
%bleach-time.
%The iterable_parameters are provided in the TARDIS GUI
iterable_parameters = parameters{1};

%% Get other parameters out, which might be important for your fitting
%Get settings information - contains settings like localisation
%uncertainty, frame-time, etc
settings_parameters = parameters{2};

%Get the current tau value (temporal shift)
tau_value = parameters{3};

%% Custom PDF area starts here
%Get readable parameters
D = iterable_parameters(1)/1e12;
loc_unc = settings_parameters.loc_unc;
strobo_frame_time = settings_parameters.strobo_frame_time*tau_value;
frame_time = settings_parameters.frame_time*tau_value;

%Following Berglund 2010
% Add loc uncertainty and motion blur coefficient to D value
Dcorr = 4*D*frame_time+2*2*loc_unc^2+4*D*(strobo_frame_time/frame_time)*(1/6)*frame_time;

% Calculate likelihood
output = xdata.*exp(-(xdata.^2)./(Dcorr));
end
