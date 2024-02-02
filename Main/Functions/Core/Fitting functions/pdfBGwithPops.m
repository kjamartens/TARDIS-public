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
function [output, extraoutput] = pdfBGwithPops(inputdata,BGcurve_interp,output_BG_alldt,pdfSettings,callfromUI,varargin)
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
%% Fix ratio of BG values for 2 pop if so chosen
% keyboard
if populations == 2
    arr_parameters = cell2mat(varargin{1}); %Transform the parameters from cell to array
    %If they are fixed, make f_pop_1 the same as the input var
    %If they are fixed, make f_pop_1 the same as the input var
    if fixRatios2pop
        if fitWithBleach
            %In case of bleach, ensure that bleach half-time is the same of
            %the populations, not the fraction at every dt
            arr_parameters = [arr_parameters arr_parameters(end)];
            f_pop_1 = arr_parameters(3:3+size(size_dt,2)-1);
        else
            for dt = 1:size(size_dt,2)
                f_pop_1(dt) = arr_parameters(3);
            end
        end
    else %if they are not fixed
        f_pop_1 = arr_parameters(3:end);
    end
end

%% Calculate the 'regular' BG values from the bleach kinetics, if used
if fitWithBleach
    arr_parameters_bleach = cell2mat(varargin{1}); %Transform the parameters from cell to array
    if populations == 1
        BGfDt = zeros(1,size(size_dt,2));
        BGfDt(1) = arr_parameters_bleach(2);
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
        arr_parameters = [arr_parameters_bleach(1) fliplr(BGfDt)];
    elseif populations == 2
        f_pop_orig = f_pop_1(1); %This is the BG ratio between the two pops
        %First calculations for the BG fractions
        BGfDt = zeros(1,size(size_dt,2));
        BGfDt(1) = arr_parameters_bleach(4);
        if fixRatios2pop
            bleachtimeframes(1) = arr_parameters_bleach(end)./frame_time;
            bleachtimeframes(2) = arr_parameters_bleach(end)./frame_time;
        else
            bleachtimeframes(1) = arr_parameters_bleach(end-1)./frame_time;
            bleachtimeframes(2) = arr_parameters_bleach(end)./frame_time;
        end
        %Now we get the decay based on the bleach time only - this still
        %needs to be corrected for nr of JDs and BG fraction at dt = 1
        BGfDt_decay(:,1) = multiExpFitTARDISBleach(bleachtimeframes(1),[1:size(size_dt,2)],[1:size(size_dt,2)]);
        BGfDt_decay(:,2) = multiExpFitTARDISBleach(bleachtimeframes(2),[1:size(size_dt,2)],[1:size(size_dt,2)]);
        BGfDt_decay = BGfDt_decay./size_dt';
        BGfDt_decay = BGfDt_decay./max(BGfDt_decay);
        %Now using this to calculate the BGfractions
        SignalFractions = (1-BGfDt(1)).*BGfDt_decay';
        %Now correct these SignalFractions for the ratio between the
        %populations
        SignalFractionSingle = SignalFractions(1,:)*(1/(1+f_pop_orig))+SignalFractions(2,:)*(f_pop_orig/(1+f_pop_orig));
        BGfDt = 1-SignalFractionSingle;

        f_pop_1 = (SignalFractions(2,:)./SignalFractions(1,:))*f_pop_orig;

        %Restore varargin
        arr_parameters = [arr_parameters_bleach(1:2) f_pop_1 fliplr(BGfDt)];
    end
else %if 'normal', no-bleach fitting is done
    arr_parameters = cell2mat(varargin{1}); %Transform the parameters from cell to array
end
%% Pull out information on fitting params from varargin
%Note: removed loc unc from varargin and use D in um2/s as input
arrsizerequired = (populations+1+populations-1+size(size_dt,2)-1); %Get the required size of the arr, based on nr of populations and size_dt
if populations == 1 %Exception for 1 population
    arrsizerequired = arrsizerequired;
end
%Check if the correct array size is entered, else give a warning
if size(arr_parameters,2) ~= arrsizerequired
%     disp('Wrong entry params - exception for 2pop')
end
%Give the parameters more tangible names. 
if populations == 1
    D1 = arr_parameters(1)/1e12;
    if verbose
        %We use a custom function to ensure that we always have the correct
        %number of digits
        string = ['Fitting params: D1: ' roundNrDigits(D1*1e12,2,2) ' um2/s | loc_unc: ' roundNrDigits(loc_unc*1e9,2,0) ' nm'];
        if fitWithBleach
            string = [string ' | Bleach half-time: ' roundNrDigits(arr_parameters_bleach(end),1,3) 's'];
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
end
if populations == 2
%     disp(['R: ' num2str(f_pop_1) ' Bl: ' num2str(bleachtimeframes)])
    D1 = arr_parameters(1)/1e12;
    D2 = arr_parameters(2)/1e12;
%     f_pop_1 = arr_parameters(3); %This is now done earlier
    if verbose
        %We use a custom function to ensure that we always have the correct
        %number of digits
        string = ['Fitting params: D1: ' roundNrDigits(D1*1e12,2,2) ' um2/s | D2: ' roundNrDigits(D2*1e12,2,2) ' um2/s | loc_unc: ' roundNrDigits(loc_unc*1e9,2,0) ' nm'];
        if fitWithBleach %This should be optimized with second bleach val
            if fixRatios2pop
                string = [string ' | Ratio_dt1: ' roundNrDigits(100-100*f_pop_1(1)/(f_pop_1(1)+1),2,1) ':' roundNrDigits(100*f_pop_1(1)/(f_pop_1(1)+1),2,1) ' | Bleach half-time 1: ' roundNrDigits(arr_parameters_bleach(end),1,3) 's' ' | Bleach half-time 2: ' roundNrDigits(arr_parameters_bleach(end),1,3) 's'];
            else
                string = [string ' | Ratio_dt1: ' roundNrDigits(100-100*f_pop_1(1)/(f_pop_1(1)+1),2,1) ':' roundNrDigits(100*f_pop_1(1)/(f_pop_1(1)+1),2,1) ' | Bleach half-time 1: ' roundNrDigits(arr_parameters_bleach(end-1),1,3) 's' ' | Bleach half-time 2: ' roundNrDigits(arr_parameters_bleach(end),1,3) 's'];
            end
         else
            if fixRatios2pop 
                string = [string ' | Ratio: ' roundNrDigits(100-100*f_pop_1(1)/(f_pop_1(1)+1),2,1) ':' roundNrDigits(100*f_pop_1(1)/(f_pop_1(1)+1),2,1)];
            else
                for dt = 1:size(size_dt,2)
                    string = [string ' | Ratio at dt ' num2str(dt) ': ' roundNrDigits(100-100*f_pop_1(dt)/(f_pop_1(dt)+1),2,1) ':' roundNrDigits(100*f_pop_1(dt)/(f_pop_1(dt)+1),2,1)];
                end
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
            else %Console output
                consoleOutputReplaceMsg(string);
            end
        catch
            fprintf([string '\n']);
        end
    end
end
if populations == 3
    D1 = arr_parameters(1)/1e12;
    D2 = arr_parameters(2)/1e12;
    D3 = arr_parameters(3)/1e12;
%     loc_unc = 15e-9;%arr_parameters(4);  %Originally, this should be dependant on arr_parameters(42)
    f_pop_1 = arr_parameters(4);
    f_pop_2 = arr_parameters(5);
    if verbose
        fprintf('%.2s - %.2s - %.2s - %.2s - %.2s - %.2s - %.2s\n',D1,D2,D3,loc_unc,f_pop_1,f_pop_2)
    end
end
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
    output_singlePop1 = pdfSinglePopFunction(inputdata(startpos(dt):endpos(dt)),D1,0,loc_unc,strobo_frame_time,frame_time*dt);
    %     output_singlePop1 = pdfSinglePopFunction(inputdata(startpos(dt):endpos(dt)),D1,0,loc_unc,1,BGcurve_interp,frame_time*dt);
    if populations > 1
        %If more than one population is used, call a second
        %pdfSinglePopFunction with relative intensity f_pop_1. Note that
        %the only other changed variable is D2 instead of D1
        %         output_singlePop2 = pdfSinglePopFunction(inputdata(startpos(dt):endpos(dt)),D2,0,loc_unc,f_pop_1,BGcurve_interp,frame_time*dt);
        output_singlePop2 = pdfSinglePopFunction(inputdata(startpos(dt):endpos(dt)),D2,0,loc_unc,strobo_frame_time,frame_time*dt);
    end
    if populations > 2
        %If more than two populations are used, call a third
        %pdfSinglePopFunction with relative intensity f_pop_3. Note that
        %the only other changed variable is D3 instead of D1
        output_singlePop3 = pdfSinglePopFunction(inputdata(startpos(dt):endpos(dt)),D3,0,loc_unc,strobo_frame_time,frame_time*dt);
        %         output_singlePop3 = pdfSinglePopFunction(inputdata(startpos(dt):endpos(dt)),D3,0,loc_unc,f_pop_2,BGcurve_interp,frame_time*dt);
    end
  
    
    %Normalize the areas of all populations (BG and up to 3 populations).
    %The normalization of the area basically scales the intensity of all
    %inputs via pseudo-integration, so that the ratio of the populations
    %and bg can be accurately set
    %Note that the function is different for normalizing area of output_BG.
    %Moreover, the '0,0' values are irrelevant
    output_BG_orig = output_BG;
    output_singlePop1_orig = output_singlePop1;
%     keyboard
    [output_BG] = normalize_area_IntFit(inputdata(startpos(dt):endpos(dt)),output_BG,'pdfBGFunction',0,0,1,BGcurve_interp,0,0);
    [output_singlePop1] = normalize_area_IntFit(inputdata(startpos(dt):endpos(dt)),output_singlePop1,'pdfSinglePopFunction',D1,loc_unc,1,BGcurve_interp,strobo_frame_time,frame_time*dt);
    
    if populations > 1
        [output_singlePop2] = normalize_area_IntFit(inputdata(startpos(dt):endpos(dt)),output_singlePop2,'pdfSinglePopFunction',D2,loc_unc,f_pop_1,BGcurve_interp,strobo_frame_time,frame_time*dt);
    end
    if populations > 2
        [output_singlePop3] = normalize_area_IntFit(inputdata(startpos(dt):endpos(dt)),output_singlePop3,'pdfSinglePopFunction',D3,loc_unc,f_pop_2,BGcurve_interp,strobo_frame_time,frame_time*dt);
    end
    %Now create the full output for this dt bin, scaling every population
    %on their respective fraction (f_BG_dt and f_pop_1 and/or f_pop_2)
    if populations == 1 %Single population only cares about f_BG_dt
        outputdt{dt} = (1*output_singlePop1).*(1-f_BG_dt{dt})+(f_BG_dt{dt}*output_BG);
    elseif populations == 2 %Two populations scales with f_BG_dt and f_pop_1
        %General formula: ([Population_ratio * population_output]_for_all_i)*Ratio_signal+[Background output]*Ratio_background
        outputdt{dt} = ((1/(1+f_pop_1(dt)))*output_singlePop1+...
            (f_pop_1(dt)/(1+f_pop_1(dt)))*output_singlePop2).*(1-f_BG_dt{dt})+...
            (f_BG_dt{dt}*output_BG);
    elseif populations == 3
        %General formula: ([Population_ratio * population_output]_for_all_i)*Ratio_signal+[Background output]*Ratio_background
        outputdt{dt} = ((1/(1+f_pop_1+f_pop_2))*output_singlePop1+...
            (f_pop_1/(1+f_pop_1+f_pop_2))*output_singlePop2+...
            (f_pop_2/(1+f_pop_1+f_pop_2))*output_singlePop3).*(1-f_BG_dt{dt})+...
            (f_BG_dt{dt}*output_BG);
    end
    
    %Remove negative and zero values
    outputdt{dt}(outputdt{dt}<=0) = 1e-9;
    outputdt{dt}(isnan(outputdt{dt})) = 1e-9;

%     disp('EXPERIMENTAL FEATURE; pdfBGwithPops Line 320~!\n')
%     if dt < 4
%         outputdt{dt} = ones(size(outputdt{dt}));
%     end
%     
    %Give an extraoutput cell with some information for further use.
    if populations == 1
        extraoutput{:,dt} = {[output_BG],[output_singlePop1]};
    elseif populations == 2
        extraoutput{:,dt} = {[output_BG],[output_singlePop1],[output_singlePop2]};
    elseif populations == 3
        extraoutput{:,dt} = {[output_BG],[output_singlePop1],[output_singlePop2],[output_singlePop3]};
    end
end %end loop over all dt's

%After all dt's, create a final array with the outputdt arrays below one
%another.
output = zeros(sum(size_dt),1);
for dt = 1:size(size_dt,2)
    output(startpos(dt):endpos(dt)) = outputdt{dt};
end
end
