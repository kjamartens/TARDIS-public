%% Background likelihood of input data - LIN variant
% Obtain the likelihood of xdata belonging to the 'background' curve.
% Basically, the interpolated BG curve is used as a lookup table for every
% entry in xdata, and a likelihood is returned.
%---------------------------------------------------------
% Required inputs
% xdata:                All input values of the JD data
% BGcurve_interp:       Interpolated background curve from
%                       interpolate_BGCurve.m
%
% Output
% output                Likelihood of every xdata to belong to the BG
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [output] = pdfBGFunction(xdata,BGcurve_interp)
%Get more exact interpolated value for xdata on the BGcurve_interp list
%For now (maybe for ever) take yvalue based on the two bins
%surrounding it, taking distance into account
%Find the localization on the BGcurve_interp xlist
dist_xstep = BGcurve_interp(2,1)-BGcurve_interp(1,1);
start_xstep = BGcurve_interp(1,1);
xloc = xdata./dist_xstep+start_xstep+1;
%Look at the bins next to it (floor/ceil of xloc) and weigh it
%based on the closeness (mod(xloc,1)).
ydata_BG_interpfound = zeros(size(xdata,1),1);%allocate space
try
    %Some float error fixing
    xloc(ceil(xloc)>size(BGcurve_interp,1)) = xloc(ceil(xloc)>size(BGcurve_interp,1))-1e-2;
    xloc(floor(xloc)<0) = xloc(floor(xloc)<0)+1e-14;
    ydata_BG_interpfound = BGcurve_interp(floor(xloc),2).*(1-mod(xloc,1))+...
        BGcurve_interp(ceil(xloc),2).*mod(xloc,1);
catch
    disp('Error pdfBGfunction')
    keyboard %Error catching
end

output = ydata_BG_interpfound;
end
