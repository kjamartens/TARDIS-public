%% Normalization of the histogram area - Signal+BG variant
% Area normalization is important for MLE
% Calculate the area by pseudo-integration. Calculate the area for
% many small bins and correct output for the sum of those bins.
% In principle, the area of the curve specified by the parameters (D,
% loc_unc, f) is calculated, and the input is normalized for the sum of this
% area
%---------------------------------------------------------
% Required inputs
% xdata:                All input values of the JD data
% input:                Likelihood from e.g. pdfSinglePopFunction_log
% functionp:            Name of function to be fitted
% D:                    Diffusion coefficient (m2/s)
% loc_unc:              Localization uncertainty in m
% f:                    Unused
% BGcurve_interp:       Interpolated background curve from
%                       interpolate_BGCurve.m
% frame_time:           Frame time in s
% strobo_frame_time     Duration that strobo-light is on. If this is
%                       set to -1, the value of the frame will be used instead.
%
% Output
% output                Likelihood of every xdata to belong to these params
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [output] = normalize_area_IntFit(xdata,input,functionp,D,loc_unc,f,BGcurve_interp,strobo_frame_time,frame_time)
%Make the bins
area_calc_bin_size = (max(xdata)-min(xdata))/1000; 
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
    mean_y_pos_bin = pdfSinglePopFunction(mean_x_pos_bin,D,0,loc_unc,strobo_frame_time,frame_time);
end
%Get total area
area_bin = mean_y_pos_bin*area_calc_bin_size;
totarea = sum(area_bin);
%Correct for total area
output = input./totarea;
% % keyboard
% % % %%
% % figure(9);clf(9);
% % plot(mean_x_pos_bin,mean_y_pos_bin,'k-x')
% % % set(gca,'XScale','log')
% % sum(input./totarea)
end