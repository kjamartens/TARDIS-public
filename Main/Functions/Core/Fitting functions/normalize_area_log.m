%% Normalization of the histogram area - LOG variant
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
% D:                    Diffusion coefficient (m2/s) OR BG interpolated
%                       curve when functionp is BG
% offset:               Offset present in data - important when
%                       BG-subtracted data is not lying at y=0
% loc_unc:              Localization uncertainty in m
% frame_time:           Frame time in s
% norm_bins:            Number of bins for normalization, normally 10k or
%                       so
% strobo_frame_time     Duration that strobo-light is on. If this is
%                       set to -1, the value of the frame will be used instead.
% functionp:            Function to be called
%
% Output
% output                Likelihood of every xdata to belong to these params
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
%(xdata,input,functionp,D,loc_unc,f,BGcurve_interp,frame_time)
function [output] = normalize_area_log(xdata,input,functionp,DorBG,loc_unc,offset,norm_bins,strobo_frame_time,frame_time)
% keyboard
% norm_bins_orig = norm_bins;
% norm_bins = 10000;
%Make the bins
separataion_calc_bin_size_log = (log10(max(xdata))-log10(min(xdata)))/norm_bins;
%Specifiy x-range it should be checked for
check_range_log = log10(min(xdata))+separataion_calc_bin_size_log/2:separataion_calc_bin_size_log:(log10(max(xdata))-separataion_calc_bin_size_log/2); 
%Get middle x position of every bin
mean_x_pos_bin_logtolin = 10.^((check_range_log+separataion_calc_bin_size_log/2)'); 

logspacing = logspace(log10(min(xdata)),log10(max(xdata)),norm_bins+1);
area_calc_bins_log = logspacing(2:end)-logspacing(1:end-1);
%Get area of all bins via the same function as in the regular probability
%calculation. This gives the total area of the function
%Different function for BG or population

%This following line works:
% mean_y_pos_bin = pdfSinglePopFunction_log(mean_x_pos_bin_logtolin,D,offset,loc_unc,frame_time);
if strcmp(functionp,'pdfBGFunction')
    mean_y_pos_bin = pdfBGFunction_log(mean_x_pos_bin_logtolin,DorBG);
elseif strcmp(functionp,'pdfSinglePopFunction')
    mean_y_pos_bin = pdfSinglePopFunction_log(mean_x_pos_bin_logtolin,DorBG,0,loc_unc,strobo_frame_time,frame_time);
%     mean_y_pos_bin = pdfSinglePopFunction(mean_x_pos_bin_logtolin,D,loc_unc,f,BGcurve_interp,frame_time);
end

%Get total area
area_bin = mean_y_pos_bin.*area_calc_bins_log';
totarea = sum(area_bin);
% % keyboard
% % %%
% % figure(9);clf(9);
% % plot(mean_x_pos_bin_logtolin,mean_y_pos_bin,'k-x')
% % set(gca,'XScale','log')
% % sum(input./totarea)


%Correct for total area
output = input./totarea;
end