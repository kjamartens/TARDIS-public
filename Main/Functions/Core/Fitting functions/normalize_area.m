%% Normalization of the histogram area - LIN variant
% Area normalization is important for MLE
% Calculate the area by pseudo-integration. Calculate the area for
% many small bins and correct output for the sum of those bins.
% In principle, the area of the curve specified by the parameters (D,
% loc_unc, f) is calculated, and the input is normalized for the sum of this
% area
%---------------------------------------------------------
% Required inputs
% xdata:                All input values of the JD data
% input:                Likelihood from e.g. pdfSinglePopFunction
% D:                    Diffusion coefficient (m2/s)
% offset:               Offset present in data - important when
%                       BG-subtracted data is not lying at y=0
% loc_unc:              Localization uncertainty in m
% frame_time:           Frame time in s
% norm_bins:            Number of bins for normalization, normally 10k or
%                       so
% functionp:            Name of the fitting function
%
% Output
% output                Likelihood of every xdata to belong to these params
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [output] = normalize_area(xdata,input,D,offset,loc_unc,frame_time,norm_bins,strobo_frame_time,functionp)
% Make the bins
area_calc_bin_size = (max(xdata)-min(xdata))/norm_bins; 
%Specifiy x-range it should be checked for
check_range = min(xdata):area_calc_bin_size:(max(xdata)-area_calc_bin_size); 
%Get middle x position of every bin
mean_x_pos_bin = (check_range+area_calc_bin_size/2)'; 
%Get area of all bins via the same function as in the regular probability
%calculation. This gives the total area of the function
%Different function for BG or population
if strcmp(functionp,'pdfBGFunction')
    mean_y_pos_bin = pdfBGFunction(mean_x_pos_bin,BGcurve_interp);
elseif strcmp(functionp,'pdfSinglePopFunction')
    mean_y_pos_bin = pdfSinglePopFunction(mean_x_pos_bin,D,offset,loc_unc,strobo_frame_time,frame_time);
end
%Get total area
area_bin = mean_y_pos_bin*area_calc_bin_size;
totarea = sum(area_bin);
%Correct for total area
output = input./totarea;
end