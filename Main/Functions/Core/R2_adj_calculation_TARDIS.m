%% Calculation of Olkin-Pratt adjusted R-squared from model fit
function [R2_adj, R2_adj_nopop] = R2_adj_calculation_TARDIS(maxdist, dt_arr,  inputarr, bgbinningnr, JDarrBG, pdfSettings, callfromUI, parameters)
% Required input:
% maxdist: in m
% dt_arr: length for every dt
% inputarr: array of JDs of signal
% bgbinningnr: nr of bins of BG
% JDarrBG: array of JDs of BG
% pdfSettings, struct containing at least:
    % frame_time: in s
    % loc_unc: in m
    % size_dt
    % populations: nr of populations
    % fitWithBleach: Boolean whether fitted with Bleach
    % fixRatios2pop: Boolean whether 2 pop bleach is fixed
% callfromUI: detailed of UI information
% parameters: output of MLE fit

% Output:
% R2_adj: OP-adjusted R2 of model fit
% R2_adj_noPop: OP-adjusted R2 of only BG fit

%Set verbose to 0 for later
pdfSettings.verbose = 0;

% Get a histogram of the inputdata
nrbins_forR2 = 300;
hist_binEdges = [linspace(0,maxdist,nrbins_forR2+1)]';
binval = zeros(nrbins_forR2*max(dt_arr),1);
hist_R2_binmids = zeros(nrbins_forR2*max(dt_arr),1);
for dt = dt_arr
    hist_R2 = histcounts(inputarr(1+sum(pdfSettings.size_dt(1:dt-1)):sum(pdfSettings.size_dt(1:dt))),hist_binEdges,'Normalization','Probability');
    hist_R2_binmids(1+(dt-1)*nrbins_forR2:(dt)*nrbins_forR2) = hist_binEdges(1:end-1)+hist_binEdges(2)/2;
    binval(1+(dt-1)*nrbins_forR2:(dt)*nrbins_forR2) = hist_R2;
end

maxdtbins = max(dt_arr);
xdataplotsingledt = [linspace(0,maxdist,nrbins_forR2)]'; %Create bins used for visualisation, based on figurebinnr
xdataplot = repmat(xdataplotsingledt,maxdtbins,1); %Extend these bins for every dtbin, and put them below each other
size_dt_plot = ones(1,maxdtbins)*(nrbins_forR2); %Make a 1-by-dt array, each holding the size of the bins, to later split up xdataplot again
pdfSettings.size_dt = size_dt_plot;
BGcurve_interp = interpolate_BGCurve(maxdist,bgbinningnr,(maxdist/bgbinningnr),JDarrBG);
output_BG_alldt = pdfBGFunction(xdataplot,BGcurve_interp); %Get BG values for all x positions (repeating for dt positions)
%Get signal values for all x positions and all dt positions
[modelData,extramodelData] = pdfBGwithPops(xdataplot,BGcurve_interp,output_BG_alldt,pdfSettings,callfromUI,{parameters});

modelData = modelData./sum(modelData)*max(dt_arr);
modelData_noPop = [];
for dt = dt_arr
    modelData_noPop_dt{dt} = extramodelData{dt}{1}./sum(extramodelData{dt}{1});
    modelData_noPop = [modelData_noPop;modelData_noPop_dt{dt}];
end
SStot = sum((binval-mean(binval)).^2);
SSres = sum((binval-modelData).^2);

R2 = 1-SSres/SStot;

R2_adj = Olkin_Pratt_Estimator(R2,max(size(hist_R2_binmids)),size(parameters,2),5);

SStot_nopop = sum((binval-mean(binval)).^2);
SSres_nopop = sum((binval-modelData_noPop).^2);

R2_nopop = 1-SSres_nopop/SStot_nopop;
R2_adj_nopop = Olkin_Pratt_Estimator(R2_nopop,max(size(hist_R2_binmids)),0,5);

%Keeping the code to show the R2 here:
% fprintf('\n')
% fprintf('R2:\t\t\t\t%.10f\n',R2);
% fprintf('R2_noPop:\t\t%.10f\n',R2_nopop);
% fprintf('R2adj:\t\t\t%.10f\n',R2_adj);
% fprintf('R2adj_noPop:\t%.10f\n',R2_adj_nopop);