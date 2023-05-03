%% Background likelihood of input data - LOG variant
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
function [output] = pdfBGFunction_log(xdata,BGcurve_interp)
%Get more exact interpolated value for xdata on the BGcurve_interp list
%For now (maybe for ever) take yvalue based on the two bins
%surrounding it, taking distance into account
%Find the localization on the BGcurve_interp xlist
startlogpos = log10(BGcurve_interp(1,1));
endlogpos = log10(BGcurve_interp(end,1));
logstep = (endlogpos-startlogpos)./(size(BGcurve_interp,1)-1);
xloclog = (log10(xdata)-startlogpos)./logstep+1;

%Look at the bins next to it (floor/ceil of xloc) and weigh it
%based on the closeness (mod(xloc,1)).
ydata_BG_interpfound = zeros(size(xdata,1),1);%allocate space
try
    %Some float error fixing
    xloclog(ceil(xloclog)>size(BGcurve_interp,1)) = xloclog(ceil(xloclog)>size(BGcurve_interp,1))-1e-14;
    xloclog(floor(xloclog)<0) = xloclog(floor(xloclog)<0)+1e-14;
    %escape 0 and too large values
    xloclog(floor(xloclog)<1) = 1;
    xloclog(ceil(xloclog)>size(BGcurve_interp,1)) = size(BGcurve_interp,1);
    ydata_BG_interpfound = (BGcurve_interp(floor(xloclog),2).*(1-mod(xloclog,1))+...
        BGcurve_interp(ceil(xloclog),2).*mod(xloclog,1))./xdata;
catch
    keyboard %Error catching
end

output = ydata_BG_interpfound;
end