%% Main population fitting method for TARDIS
% Supposed to fit data with populations and background (intra-track
% linkages).
%---------------------------------------------------------
% Required inputs
% inputdata          The JD data that should be fitted/compared
% BGcurve_interp     L-x-2 matrix with 1st column x-positions, 2nd column
%                    y-positions of the BGcurve. Also see interpolate_BGCurve, where this is an
%                    output of. Normally run outside this function to speed up
%                    analysis time, but interpolate_BGCurve could be inside
%                    this function
% output_BG_alldt    Non-normalized probability of all inputdata values
%                    belonging to the BG. This is normally created via the
%                    function pdfBGFunction outside this script.
% pdfSettings        Struct containing the following (at least):
%      frame_time         Frame time in seconds
%      populations        Number of populations
%      size_dt            1-by-dt array with number of values of inputdata that
%                         correspond to each dt. Therefore, sum(size_dt) should
%                         equal the size of inputdata.
%      verbose            True for verbose output, False for silent operation
%      fitWithBleach      Boolean: true if the BG/populations are fitted with
%                         bleach curves, false if the BG/populations are fitted
%                         with 'normal' arbitrary values
%      fixRatios2pop      Boolean wether or not Ratios of 2 pops should be
%                         fixed over all DTs
%      strobo_frame_time  Duration that strobo-light is on. If this is
%                         set to -1, the value of the frame will be used instead.
% callfromUI         Information about the UI, or empty if run from
%                       console.
% varargin           {1,1}-cell containing the parameters for the model. See
%                    earlier remarks for size and contents of parameters at
%                    the start of the wrapper function.
%
% Obtained outputs:
% output            Array of size inputdata with likelihood of it belonging
%                   to the parameters specified in varargin
% extraoutput       Information required for visualistion
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [output, extraoutput] = CustomPDFwithBG(inputdata,BGcurve_interp,output_BG_alldt,pdfSettings,callfromUI,PDFname,varargin)
%% Extract all variables from pdfSettings
frame_time = pdfSettings.frame_time;
loc_unc = pdfSettings.loc_unc;
populations = pdfSettings.populations;
size_dt = pdfSettings.size_dt;
verbose = pdfSettings.verbose;
fitWithBleach = pdfSettings.fitWithBleach;
fixRatios2pop = pdfSettings.fixRatios2pop;
freefit_locunc = pdfSettings.freefit_locunc;
strobo_frame_time = pdfSettings.strobo_frame_time;
populations = 1;
% fitWithBleach = 1;

%% Calculate the 'regular' BG values from the bleach kinetics, if used
if fitWithBleach
    arr_parameters_bleach = cell2mat(varargin{1}); %Transform the parameters from cell to array
    
    BGfDt = zeros(1,size(size_dt,2));
    BGfDt(1) = arr_parameters_bleach(end-1);
    bleachtimeframes = arr_parameters_bleach(end)./frame_time;
    %Now we get the decay based on the bleach time only - this still
    %needs to be corrected for nr of JDs and BG fraction at dt = 1
    %*(max(size_dt)/min(size_dt))^2
    BGfDt_decay = multiExpFitTARDISBleach(bleachtimeframes,[1:size(size_dt,2)],[1:size(size_dt,2)]);
    BGfDt_decay = BGfDt_decay./size_dt';
    BGfDt_decay = BGfDt_decay./max(BGfDt_decay);
    %Now using this to calculate the BGfractions
    SignalFractions = (1-BGfDt(1)).*BGfDt_decay';
    BGfDt = 1-SignalFractions;
    
    %Restore varargin
    arr_parameters = [arr_parameters_bleach(1:end-2) fliplr(BGfDt)];
else %if 'normal', no-bleach fitting is done
    arr_parameters = cell2mat(varargin{1}); %Transform the parameters from cell to array
end
%% Pull out information on fitting params from varargin
%f_BG_dt is fraction of BG for every dt
for k = 1:size(size_dt,2)
    f_BG_dt{k} = arr_parameters(end-(k-1));
end
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
%% Calculate probabilities that every input value corresponds to the model function
%Calculate the probability of all entries in inputdata that it belongs to
%the model described by the parameters set up above, and the
%BGcurve_interp. The fraction of BGcurve_interp is described via f_BG_dt.

for dt = 1:size(size_dt,2) %Loop over all dt inputs
    %Select output from output_BG_alldt that correspond to this dt. This is
    %basically the likelihood that a point belongs to the background. This
    %is later scaled according to f_BG_dt
    output_BG = output_BG_alldt(startpos(dt):endpos(dt));
    %Get the likelihood that the input points belong to either of the
    %populations (pop1 or pop2, maybe a pop3 should be added later)
    %The used ratio of populations is 1:f_pop_1:f_pop_2 if 2 populations
    %are used.
    
    %When at least one population is used, call pdfSinglePopFunction with
    %relative intensity 1

    %Create the parameters for this run
    %Last 2 entries of varargin are BG and bleach settings, respectively.
    params = {arr_parameters(1:end-size(size_dt,2)),pdfSettings,dt};

    %Evaluate
    output_singlePop1 = eval([PDFname '(inputdata(startpos(dt):endpos(dt)),params)']);%D1,0,loc_unc,strobo_frame_time,frame_time*dt);
    
    %Normalize the areas of all populations (BG and up to 3 populations).
    %The normalization of the area basically scales the intensity of all
    %inputs via pseudo-integration, so that the ratio of the populations
    %and bg can be accurately set
    %Note that the function is different for normalizing area of output_BG.
    %Moreover, the '0,0' values are irrelevant
    output_BG_orig = output_BG;
    output_singlePop1_orig = output_singlePop1;

    %Noramlise - same PDF and parameters
    [output_BG] = normalize_area_IntFit(inputdata(startpos(dt):endpos(dt)),output_BG,'pdfBGFunction',0,0,1,BGcurve_interp,0,0);
    [output_singlePop1] = normalize_area_IntFit(inputdata(startpos(dt):endpos(dt)),output_singlePop1,PDFname,0,0,params,BGcurve_interp,0,0);
  
    %Now create the full output for this dt bin, scaling every population
    %on their respective fraction (f_BG_dt and f_pop_1 and/or f_pop_2)
    outputdt{dt} = (1*output_singlePop1).*(1-f_BG_dt{dt})+(f_BG_dt{dt}*output_BG);
   
    
    %Remove negative and zero values
    outputdt{dt}(outputdt{dt}<=0) = 1e-9;
    outputdt{dt}(isnan(outputdt{dt})) = 1e-9;

    %Give an extraoutput cell with some information for further use.
    extraoutput{:,dt} = {[output_BG],[output_singlePop1]};
end %end loop over all dt's
%% Verbose output
%Give the parameters more tangible names. 
if verbose
    %We use a custom function to ensure that we always have the correct
    %number of digits
    string = ['CustomPDF - ' PDFname ' - Fitting params: '];
    for i = 1:size(arr_parameters(1:end-size(size_dt,2)),2)
        string = [string roundNrDigits(arr_parameters(i),1,3) ' | '];
    end
    if fitWithBleach
        string = [string 'Bleach half-time: ' roundNrDigits(arr_parameters_bleach(end),1,3) 's'];
    else
        string = [string 'intra-emitter fractions: '];
        for i = size(arr_parameters,2)-size(size_dt,2)+1:size(arr_parameters,2)
            string = [string roundNrDigits(arr_parameters(i),1,3) ' | '];
        end
    end
    %Output two or three dots to show it's ongoing
    if rand() < 0.5
        string = [string '..'];
    else
        string = [string '...'];
    end
    try
        if ~(isempty(callfromUI))
            %UI output
            %Testing writing somewhere
            %Check if last line needs to be removed
            if contains(char(callfromUI.TextArea.Value(1)),'Fitting params:')
                replaceLastTextLine(string,callfromUI);
            else
                dispUIorCommandWindow(string,callfromUI);
            end
        else %Console output - prevmsglen is retained in global var
            consoleOutputReplaceMsg(string);
        end
    catch
        fprintf([string '\n']);
    end
end
%% Final output
%After all dt's, create a final array with the outputdt arrays below one
%another.
output = zeros(sum(size_dt),1);
for dt = 1:size(size_dt,2)
    output(startpos(dt):endpos(dt)) = outputdt{dt};
end

end
