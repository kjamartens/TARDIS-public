%% Custom PDF fitting routine - Single diffusive species under rod-shape confinement
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
function output = CustomPDF_CellShapeConfinement_1pop(xdata,parameters)
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

%Correct further for confinement - based on anaDDA
%Core idea: take ratio of cell length and cell width (on a volume-basis),
%and account for confinement in a sphere (cell caps) or circle (cell
%cylinder) based on this ratio

%We load some data about the zeros of the BesselFunction
load('ZerosOfConfinement.mat','zerovaluesCircular','zerovaluesSpherical');

%Determine the sphere radius from iterable parameters
cellRadius = iterable_parameters(2)/1e6;
cellLength = iterable_parameters(3)/1e6+2*cellRadius;

%From anaDDA
Volumecylinder = pi().*cellRadius.^2.*(cellLength-2.*cellRadius);
Volumesphere = 4./3*pi().*cellRadius.^3;
RatioCircleConfinement = Volumecylinder/(Volumecylinder+Volumesphere);
RatioSphereConfinement = Volumesphere/(Volumecylinder+Volumesphere);

%Number of zeros to assess
nritsum = 10000;

%Calculate the sum-part of the Dcorr factor Spherical
p1 = exp(-(zerovaluesSpherical(1:nritsum).^2.*frame_time.*Dcorr)./(cellRadius^2));
p2 = (1./(zerovaluesSpherical(1:nritsum).^2.*(zerovaluesSpherical(1:nritsum).^2-2)));
sumSpherical = sum(p1.*p2);

%Calculate the sum-part of the Dcorr factor Circular
p1 = exp(-(zerovaluesCircular(1:nritsum).^2.*frame_time.*Dcorr)./(cellRadius^2));
p2 = (1./(zerovaluesCircular(1:nritsum).^2.*(zerovaluesCircular(1:nritsum).^2-1)));
sumCircular = sum(p1.*p2);

%Get Dcorr in x and y, confined for spherical and cylindrical parts
% Assumption is that cell's long axis is x
Dcorr_x_spherical = (cellRadius^2/(6*frame_time))*(6/5-12*sumSpherical); %um2/s
Dcorr_x_cylinder = Dcorr;
Dcorr_y_spherical = (cellRadius^2/(6*frame_time))*(6/5-12*sumSpherical); %um2/s
Dcorr_y_cylinder = (cellRadius^2/(4*frame_time))*(1-8*sumCircular);
%Combine to 'observable' Dcorr in x and y, based on their ratios
Dcorr_obs_x = RatioSphereConfinement*Dcorr_x_spherical+RatioCircleConfinement*Dcorr_x_cylinder;
Dcorr_obs_y = RatioSphereConfinement*Dcorr_y_spherical+RatioCircleConfinement*Dcorr_y_cylinder;

%Comebine Dx,Dy
Dcorr_cellShape = sqrt(Dcorr_obs_x*Dcorr_obs_y);

% Calculate likelihood
output = xdata.*exp(-(xdata.^2)./(4*Dcorr_cellShape*frame_time));
end
