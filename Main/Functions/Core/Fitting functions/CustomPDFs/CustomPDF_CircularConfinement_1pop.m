%% Custom PDF fitting routine - Single diffusive species under circular confinement
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
function output = CustomPDF_CircularConfinement_1pop(xdata,parameters)
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
loc_unc = settings_parameters.loc_unc;
strobo_frame_time = settings_parameters.strobo_frame_time*tau_value;
frame_time = settings_parameters.frame_time*tau_value;

%Following Berglund 2010
% Add loc uncertainty and motion blur coefficient to D value
Dcorr = (4*D*frame_time+2*2*loc_unc^2+4*D*(strobo_frame_time/frame_time)*(1/6)*frame_time)/(4*frame_time);

%Correct further for confinement
%We load some data about the zeros of the BesselFunction
load('ZerosOfConfinement.mat','zerovaluesCircular','zerovaluesSpherical');

%Determine the sphere radius from iterable parameters
sphereRadius = iterable_parameters(2)/1e6;
%Number of zeros to assess
nritsum = 10000;
%Calculate the sum-part of the Dcorr factor
p1 = exp(-(zerovaluesSpherical(1:nritsum).^2.*frame_time.*Dcorr)./(sphereRadius^2));
p2 = (1./(zerovaluesSpherical(1:nritsum).^2.*(zerovaluesSpherical(1:nritsum).^2-2)));
sum2 = sum(p1.*p2);
%Correct for circular confinement
Dcorr_spherical = (sphereRadius^2/(6*frame_time))*(6/5-12*sum2); %um2/s

% Calculate likelihood
output = xdata.*exp(-(xdata.^2)./(4*Dcorr_spherical*frame_time));
end
