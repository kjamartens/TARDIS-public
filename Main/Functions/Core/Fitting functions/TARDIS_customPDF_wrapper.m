%Wrapper function to go from TARDIS info to found parameters - should be
%changed for every custom function.

function [parameters, parametersCI] = TARDIS_customPDF_wrapper(customPDF, inputarr, bgratios, linorlogBGsubtract, BGcurve_interp,output_BG_alldt,pdfSettings,customPDF_start, customPDF_lower, customPDF_upper,callfromUI)
    %% Custom PDF routine
    %Get the correct start, lower, upper parameters based on nr of
    %populations and estimated parameters etc
    startparam = [customPDF_start fliplr(bgratios)];
    lowerbound = [customPDF_lower max(0,fliplr(bgratios)*.8)];
    upperbound = [customPDF_upper min(1,fliplr(bgratios)*1.2)];

    %Change parameters and lower/upper bound in case we fit with bleach
    %kinetics rather than fully random BG values
    [startparam, lowerbound, upperbound] = changeFitParams_populations_Bleach_customPDF(pdfSettings.fitWithBleach, size(customPDF_start,2), startparam, lowerbound, upperbound);

    %Create a customPDF
    %The used function can be found in the subscript defined here
    %(pdfBGwithPops), which has a certain set inputs, and the fitting variables
    %are in the varargin
    
    if strcmp(linorlogBGsubtract,'lin')
        pdf_Custom = @(xdata,varargin)CustomPDFwithBG(xdata,BGcurve_interp,output_BG_alldt,pdfSettings,callfromUI,customPDF,varargin);
    else
        disp('Use LIN BG subtraction!')
    end
    %perform fit
    mleoptions.Display='iter';
    %Make it harder to fit, i.e. more precise
    mleoptions.TolFun = 1e-12;
    mleoptions.TolX = 1e-12;
    mleoptions.MaxFunEvals = 10000;
    mleoptions.MaxIter = 3000;
    
    [parameters, parametersCI] = mle(inputarr, 'pdf',pdf_Custom,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound,'Options',mleoptions);
end