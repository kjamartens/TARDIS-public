%% Main aDDA fitting method for TARDIS
% Supposed to fit data with populations and background (intra-track
% linkages).
%---------------------------------------------------------
% Required inputs
% inputdata         The JD data that should be fitted/compared
% BGcurve_interp    L-x-2 matrix with 1st column x-positions, 2nd column
%                   y-positions of the BGcurve. Also see interpolate_BGCurve, where this is an
%                   output of. Normally run outside this function to speed up
%                   analysis time, but interpolate_BGCurve could be inside
%                   this function
% inputSettings     aDDA input: Settings structure
% rangeD            D values calculated from JD at every dt. JD is static
%                   for all dts, but D scales with JD
% fitspecies        aDDA input: nr of species to fit
% fixedspecies      aDDA input: nr of fixed species
% Dfixed            aDDA input: D-information of fixed species
% fixedparameters   aDDA input: array of size 4x3 with 1s where params are
%                   fixed, -1 elsewhere
% indexfittingparameters aDDA input: array of size 4x3 with 1s where params
%                   should be changed/fitted
% fx:               input for fittinganaDDA
% fy:               same as fx
% maxDindtracking:  input for fittinganaDDA
% frequency:        size_dt info, but matricized for multiple fits
% size_dt:          Nr JDs for every dt
% locerrorpdfcorrected_singledtbin: Pre-calculated PDF for the loc-error
% output_BG_alldt:  Likelihood of every inputdata belonging to BG
% fitWithBleach:    1 if fitting with bleach kinetics, 0 otherwise
% verbose           True for verbose output, False for silent operation
% callfromUI        Information about the UI, or empty if run from
%                       console.
% varargin          {1,1}-cell containing the parameters for the model. See
%                   earlier remarks for size and contents of parameters at
%                   the start of the wrapper function.

% Obtained outputs:
% output            Array of size inputdata with likelihood of it belonging
%                   to the parameters specified in varargin
% extraoutput       Information required for visualistion
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [output,extraoutput] = pdfAnaDDAMLE_multidt(inputdata,BGcurve_interp, Numberofframes,Frametimelist,inputsettings,rangeD, fitspecies, fixedspecies, Dfixed,fixedparameters,indexfittingparameters,fx,fy,maxDindtracking,frequency,size_dt,locerrorpdfcorrected_singledtbin,output_BG_alldt,bgbinningnr,fitWithBleach,verbose,callfromUI,varargin)
%Output of this part:
%parameters: [4x3] matrix with [0 k1 k2 D] as first line
%f_BG_dt: {1xDT} cell of BG fractions
parameters = fixedparameters;
if fitWithBleach
    try
        arr_parameters_bleach = cell2mat(varargin{1}); %Transform the parameters from cell to array
    catch
        arr_parameters_bleach = varargin{1}; %Transform the parameters from cell to array
    end
    parameters(indexfittingparameters) = arr_parameters_bleach(1:3);

    BGfDt = zeros(1,size(size_dt,2));
    BGfDt(1) = arr_parameters_bleach(4); %First 3 are aDDA input, then first BG ratio
    bleachtimeframes = arr_parameters_bleach(end)./inputsettings.frametime;
    %Now we get the decay based on the bleach time only - this still
    %needs to be corrected for nr of JDs and BG fraction at dt = 1
    BGfDt_decay = multiExpFitTARDISBleach(bleachtimeframes,[1:size(size_dt,2)],[1:size(size_dt,2)]);
    %First corrected for nr of JDs
    BGfDt_decay = BGfDt_decay./size_dt';
    BGfDt_decay = BGfDt_decay./max(BGfDt_decay);
    %Now using this to calculate the BGfractions
    SignalFractions = (1-BGfDt(1)).*BGfDt_decay';
    BGfDt = 1-SignalFractions;
    for k = 1:size(size_dt,2)
        f_BG_dt{k} = BGfDt(k);
    end
else %if 'normal', no-bleach fitting is done
    try
        parameters(indexfittingparameters) = cell2mat(varargin{1}(1:end-size(size_dt,2)));
    catch
        t = cell2mat(varargin);
        parameters(indexfittingparameters) = (t(1:end-size(size_dt,2)));
        clear t
    end
    try
        for k = 1:size(size_dt,2)
            f_BG_dt{k} = cell2mat(varargin{1}(end-k+1));
        end
    catch
        t = cell2mat(varargin);
        for k = 1:size(size_dt,2)
            f_BG_dt{k} = (t(end-k+1));
        end
        clear t
    end
end
%%
% x: Values at which to calculate the likelihood
% BGcurve_interp:    Interpolated background curve at common X-intervals
    
%Verbose text output on fitting
try
    string = '';
    string = [string 'Fitting params: k1: ' roundNrDigits(cell2mat(varargin{1}(1)),2,2) ' | k2: ' roundNrDigits(cell2mat(varargin{1}(2)),2,2) ' | D: ' roundNrDigits(cell2mat(varargin{1}(3)),2,2) ];    
    if fitWithBleach
        string = [string ' | Bleach half-time: ' roundNrDigits(cell2mat(varargin{1}(end)),1,3) 's'];
    else
        for i = 1:(size(varargin{1},2)-3)
            string = [string ' | fBG' num2str(i) ': ' roundNrDigits(cell2mat(varargin{1}(end-(i-1))),1,2)];
        end
    end
    if verbose
        if ~(isempty(callfromUI))
            %UI output
            %Testing writing somewhere
            %Check if last line needs to be removed
            if contains(char(callfromUI.TextArea.Value(1)),'Fitting params:')
                replaceLastTextLine(string,callfromUI);
            else
                dispUIorCommandWindow(string,callfromUI);
            end
        else %Console output
            consoleOutputReplaceMsg(string);
        end
    end
catch
    if verbose
        disp('...')
    end
end
%% Get parameters in a more readable form
c = parameters(:,1);
koff = parameters(:,2);
kon = parameters(:,3);
Dfree = parameters(:,4);
if inputsettings.numberofspecies<3
    c(inputsettings.numberofspecies+1:end) = 0;
end
c(1) = 1 - c(2) - c(3);

%% Separate input data on dt-step information
% For all dt possibilities, determine start and end pos based on size_dt
% array. These parameters are later used to split inputdata and
% BGcurve_interp
startpos = ones(size(size_dt,2),1); %pre-allocate
endpos = ones(size(size_dt,2),1); %pre-allocate
for k = 1:size(size_dt,2)
    if k > 1
        startpos(k) = size_dt(k-1)+startpos(k-1); %note that starpos(1) remains at 1, which should always be the case
    end
    endpos(k) = startpos(k)+size_dt(k)-1;
end

%% If one of the kinetic parameters is not fixed it can be directly calculated from the meanD
% Dmean2 = c(2)* Dfree(2)*(1-kon(2)/(koff(2)+kon(2)));
% Dmean3 = c(3)* Dfree(3)*(1-kon(3)/(koff(3)+kon(3)));

%% Used to be able to use the PDA/Stracy distribution together with MLE in MATLAB. Returns a probability for each point in the distribution.
    for dt = 1:size(size_dt,2) %Loop over all dt inputs
        try
        dataind = 1;
        %Prepare empty value
        output_singledt_anaDDA{dt} = zeros(size_dt(dt),1);
        %Calculate the aDDA combined PDF - i.e. the PDF that should fully
        %describe the data if explained by aDDA
        maxindex = numel(rangeD(dt,:));
        inputsettings.frametime = inputsettings.frametimerange(1).*dt;
        locerrorpdfcorrected = locerrorpdfcorrected_singledtbin{dt};
        %Here actual aDDA likelihood is being calculated
        [combinedpdf] = fittinganaDDA(fitspecies, koff, kon, Dfree, rangeD(dt,:), locerrorpdfcorrected, maxindex, fx, fy, maxDindtracking{dt},inputsettings,dt,c,fixedspecies,Dfixed);
        %Normalization of aDDA likelihood
        combinedpdf = combinedpdf./sum(combinedpdf);
        numberofdata = frequency(dt,fitspecies);
        ind = fitspecies*(inputdata(startpos(dt):endpos(dt))-rangeD(dt,1))/(rangeD(dt,2)-rangeD(dt,1))+1;
        outputtemp = interp1(combinedpdf(1:min(round(max(ind))+100,size(combinedpdf,1)),fitspecies),ind,'spline');%Change by KM
        output_singledt_anaDDA{dt}(dataind:dataind+numberofdata-1) = outputtemp;
        dataind = dataind+numberofdata;
        %And here we correct the likelihood based on the value of x
        output_singledt_anaDDA{dt} = output_singledt_anaDDA{dt}.*inputdata(startpos(dt):endpos(dt));

        %Get likelihood of all xpoints belonging to BG, similar to how we
        %do for aDDA
        output_BG{dt} = output_BG_alldt{dt}(startpos(dt):endpos(dt)).*inputdata(startpos(dt):endpos(dt));
        %[output] = normalize_area_log(xdata,input,functionp,DorBG,loc_unc,offset,norm_bins,frame_time)
%         [output_BG{dt}] = normalize_area_log(inputdata(startpos(dt):endpos(dt)),output_BG{dt},BGcurve_interp{dt},0,0,0,1000,'pdfBGFunction');
% keyboard

        [output_BG{dt}] = normalize_area_log(inputdata(startpos(dt):endpos(dt)),output_BG{dt},'pdfBGFunction',BGcurve_interp{dt},0,0,50000,0); %changed 1k to 50k
        %(xdata,input,functionp,DorBG,loc_unc,offset,norm_bins,strobo_frame_time,frame_time)

        %Correction factor for the BG value based on sampling factor
        maxD = (inputsettings.trackingwindow*inputsettings.pixelsize)^2/(4*inputsettings.frametimerange(1));
        k = (maxD./(rangeD(dt,2)-rangeD(dt,1)));
        c = k/maxD;
        corr_factor = (c/sum(combinedpdf)); %Correcting for precision - combinedpdf should be 1, but is not always
        output_BG{dt} = output_BG{dt}./corr_factor;
        
        %Finally, we create the total likelhood based on BG ratios
        outputdt{dt} = (1*output_singledt_anaDDA{dt}).*(1-f_BG_dt{dt})+(f_BG_dt{dt}*output_BG{dt});%output_singledt_anaDDA;
%         keyboard
        
        %Give an extraoutput cell with some information for further use.
        extraoutput{:,dt} = {[output_BG],[output_singledt_anaDDA]};
        catch
            keyboard
        end
        
    end
%After all dt's, create a final array with the outputdt arrays below one
%another.
output = zeros(sum(size_dt),1);
for dt = 1:size(size_dt,2)
    output(startpos(dt):endpos(dt)) = outputdt{dt};
end

%Limit output to be positive
output = max(1e-99,output);
end