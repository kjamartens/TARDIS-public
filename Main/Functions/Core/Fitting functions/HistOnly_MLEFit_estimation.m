%% Fitting of BG-subtracted data (HO fitting, estimation fit)
% Fit the BG-subtracted, or 'pure signal' JD histograms with the wanted
% analytical formula, to get a nice fit. Normally used as estimator for fit
% including the BG values, but could be used stand-alone.
%---------------------------------------------------------
% Required inputs
% xdata:                All input values of the JD data
% offsetpartial:        Absolute value of the baseline that is added for
%                       every dt - obtain from createReconstitutedJDarray.m
% populations:          Number of populations - 1 or 2, or 0 for aDDA
% loc_unc:              Localization uncertainty in m
% frame_time:           Frame time in s
% dtbinsize:            Number of linkages in every dt bin. Sum of this
%                       should add up to total xdata length
% norm_bins:            Number of bins used for normalization. Normally ~
%                       10k
% linorlogBGsubtract:   'lin' or 'log', depending on a linear or
%                       logarithmic BG subtract method
% debug:                Boolean to get some extra debug info
% minlogpoint:          Minimum x value in case log value is used
% verbose:              Boolean to show output during fitting
% callfromUI:           UI infromation or empty array if displayed generally
% varargin:             Actual fitting parameters. Number and meaning of
%                       parameters dependant on populations. Normally [D1]
%                       for 1 pop, [D1 D2 ratio] for 2 pop, [kon koff
%                       Dfree] for aDDA
%
% Output
% output:               Output likelihood for every xdata input
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [output] = HistOnly_MLEFit_estimation(xdata,offsetpartial,pdfSettings,dtbinsize,norm_bins,linorlogBGsubtract,debug,minlogpoint,callfromUI,varargin)
%% Extract all variables from pdfSettings
frame_time = pdfSettings.frame_time;
loc_unc = pdfSettings.loc_unc;
populations = pdfSettings.populations;
verbose = pdfSettings.verbose;
fitWithBleach = pdfSettings.fitWithBleach;
fixRatios2pop = pdfSettings.fixRatios2pop;
freefit_locunc = pdfSettings.freefit_locunc;
strobo_frame_time = pdfSettings.strobo_frame_time;
%%
arr_parameters = cell2mat(varargin{1}); %Transform the parameters from cell to array
if populations == 1
    D1 = arr_parameters(1)/1e12; %One population - only a single D is wanted
end
if populations == 2
    D1 = arr_parameters(1)/1e12;
    D2 = arr_parameters(2)/1e12;
    f_pop_1 = arr_parameters(3);
end
if populations == 3
    D1 = arr_parameters(1)/1e12;
    D2 = arr_parameters(2)/1e12;
    D3 = arr_parameters(3)/1e12;
    f_pop_1 = arr_parameters(4);
    f_pop_2 = arr_parameters(5);
end
%assuming that f.e. dtbins = [1 2 3], meaning 1,2,3 frame delay, and all
%the same x-range/xdata
%output will be [ydatadt1 ydatadt2 ydatadt3]

%Pre-allocate output
output = zeros(size(xdata));
%Loop over dt
for dt = 1:size(dtbinsize,2)
    %Adept frame time to reflect dt
    frame_time_dt = frame_time*dt;
    %calculate offset on the fly. Unsure yet about exact formula
    %     if linorlogBGsubtract == 1
    %         offset{dt} = (D1+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
    %     elseif linorlogBGsubtract == 2 %ABSOLUTELY NO CLUE WHY
    %         offset{dt} = (D1+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
    %     end
    %Get startpoint and endpoint to get the data corresponding to this dt
    %bin
    if dt > 1
        startpointxdata = sum(dtbinsize(1:dt-1))+1;
    else
        startpointxdata = 1;
    end
    endpointxdata = startpointxdata+dtbinsize(dt)-1;
    %Get the data belonging to this dt bin
    xdatadt = xdata(startpointxdata:endpointxdata);
    %Perform fitting on this dt-separated data
    if linorlogBGsubtract == 1 %Linear
        if populations > 0
            offset{dt} = (D1+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
            outputsingledt1{dt} = pdfSinglePopFunction(xdatadt,D1,offset{dt},loc_unc,strobo_frame_time,frame_time_dt);
            %Normalize for area
            outputsingledt1{dt} = normalize_area(xdatadt,outputsingledt1{dt},D1,offset{dt},loc_unc,frame_time_dt,norm_bins,strobo_frame_time,'pdfSinglePopFunction');
            %Avoid some errors (never called in good fittings)
            outputsingledt1{dt}(outputsingledt1{dt}<=0) = 1e-150;
        end
        if populations > 1
            offset{dt} = (D2+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
            outputsingledt2{dt} = pdfSinglePopFunction(xdatadt,D2,offset{dt},loc_unc,strobo_frame_time,frame_time_dt);
            %Normalize for area
            outputsingledt2{dt} = normalize_area(xdatadt,outputsingledt2{dt},D2,offset{dt},loc_unc,frame_time_dt,norm_bins,strobo_frame_time,'pdfSinglePopFunction');
            %Avoid some errors (never called in good fittings)
            outputsingledt2{dt}(outputsingledt2{dt}<=0) = 1e-150;
        end
        if populations > 2
            offset{dt} = (D3+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
            outputsingledt3{dt} = pdfSinglePopFunction(xdatadt,D3,offset{dt},loc_unc,strobo_frame_time,frame_time_dt);
            %Normalize for area
            outputsingledt3{dt} = normalize_area(xdatadt,outputsingledt3{dt},D2,offset{dt},loc_unc,strobo_frame_time,norm_bins,frame_time_dt,'pdfSinglePopFunction');
            %Avoid some errors (never called in good fittings)
            outputsingledt3{dt}(outputsingledt3{dt}<=0) = 1e-150;
        end
    elseif linorlogBGsubtract == 2 %Logarithmic
        if populations > 0
            offset{dt} = (D1+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
            outputsingledt1{dt} = pdfSinglePopFunction_log(xdatadt,D1,offset{dt},loc_unc,strobo_frame_time,frame_time_dt);
            %Normalize for area
%             keyboard
%             normalize_area_log(inputdata(startpos(dt):endpos(dt)),output_BG,'pdfBGFunction',BGcurve_interp,loc_unc,1,10000,frame_time*dt);
            %(xdata,input,functionp,DorBG,loc_unc,offset,norm_bins,frame_time)
%             outputsingledt1{dt} = normalize_area_log(xdatadt,outputsingledt1{dt},D1,offset{dt},loc_unc,frame_time_dt,norm_bins,'pdfSinglePopFunction');
            outputsingledt1{dt} = normalize_area_log(xdatadt,outputsingledt1{dt},'pdfSinglePopFunction',D1,loc_unc,offset{dt},norm_bins,strobo_frame_time,frame_time_dt);
            %Avoid some errors (never called in good fittings)
            outputsingledt1{dt}(outputsingledt1{dt}<=0) = 1e-150;
        end
        if populations > 1
            offset{dt} = (D2+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
            outputsingledt2{dt} = pdfSinglePopFunction_log(xdatadt,D2,offset{dt},loc_unc,strobo_frame_time,frame_time_dt);
            %Normalize for area
%             outputsingledt2{dt} = normalize_area_log(xdatadt,outputsingledt2{dt},D2,offset{dt},loc_unc,frame_time_dt,norm_bins,'pdfSinglePopFunction');
            outputsingledt2{dt} = normalize_area_log(xdatadt,outputsingledt1{dt},'pdfSinglePopFunction',D2,loc_unc,offset{dt},norm_bins,strobo_frame_time,frame_time_dt);
            
            %Avoid some errors (never called in good fittings)
            outputsingledt2{dt}(outputsingledt2{dt}<=0) = 1e-150;
        end
        if populations > 2
            offset{dt} = (D3+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
            outputsingledt3{dt} = pdfSinglePopFunction_log(xdatadt,D3,offset{dt},loc_unc,strobo_frame_time,frame_time_dt);
            %Normalize for area
%             outputsingledt3{dt} = normalize_area_log(xdatadt,outputsingledt3{dt},D3,offset{dt},loc_unc,frame_time_dt,norm_bins,'pdfSinglePopFunction');
            outputsingledt3{dt} = normalize_area_log(xdatadt,outputsingledt1{dt},'pdfSinglePopFunction',D3,loc_unc,offset{dt},norm_bins,strobo_frame_time,frame_time_dt);
            
            %Avoid some errors (never called in good fittings)
            outputsingledt3{dt}(outputsingledt3{dt}<=0) = 1e-150;
        end
    end
    
    %Now create the full output for this dt bin, scaling every population
    %on their respective fraction (f_BG_dt and f_pop_1 and/or f_pop_2)
    if populations == 1 %Single population only
        outputdt{dt} = outputsingledt1{dt};
    elseif populations == 2 %Two populations scales with f_BG_dt and f_pop_1
        %General formula: ([Population_ratio * population_output]_for_all_i)*Ratio_signal+[Background output]*Ratio_background
        outputdt{dt} = (1/(1+f_pop_1))*outputsingledt1{dt}+...
            (f_pop_1/(1+f_pop_1))*outputsingledt2{dt};
    elseif populations == 3
        %General formula: ([Population_ratio * population_output]_for_all_i)*Ratio_signal+[Background output]*Ratio_background
        outputdt{dt} = (1/(1+f_pop_1+f_pop_2))*outputsingledt1{dt}+...
            (f_pop_1/(1+f_pop_1+f_pop_2))*outputsingledt2{dt}+...
            (f_pop_2/(1+f_pop_1+f_pop_2))*outputsingledt3{dt};
    end
    
    %Text output
    %Give the parameters more tangible names. 
    if populations == 1
        D1 = arr_parameters(1)/1e12;
        if verbose
            %We use a custom function to ensure that we always have the correct
            %number of digits
            string = ['Fitting params: D1: ' roundNrDigits(D1*1e12,2,2) ' um2/s | loc_unc: ' roundNrDigits(loc_unc*1e9,2,0) ' nm'];
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
    if populations == 2
        D1 = arr_parameters(1)/1e12;
        D2 = arr_parameters(2)/1e12;
        f_pop_1 = arr_parameters(3);
        if verbose
            %We use a custom function to ensure that we always have the correct
            %number of digits
            string = ['Fitting params: D1: ' roundNrDigits(D1*1e12,2,2) ' um2/s | D2: ' roundNrDigits(D2*1e12,2,2) ' um2/s | loc_unc: ' roundNrDigits(loc_unc*1e9,2,0) ' nm | Ratio: ' roundNrDigits(100-100*f_pop_1/(f_pop_1+1),2,1) ':' roundNrDigits(100*f_pop_1/(f_pop_1+1),2,1)];
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

    output(startpointxdata:endpointxdata) = outputdt{dt};
end
%%
%DEBUG drawing - can later be expanded/cleaned in some nice stuff
if debug
    if populations == 1
        figure(92);clf(92);
        for dt = 1:size(dtbinsize,2)
            frame_time_dt = frame_time*dt;
            subplot(1,size(dtbinsize,2),dt)
            if dt > 1
                startpointxdata = sum(dtbinsize(1:dt-1))+1;
            else
                startpointxdata = 1;
            end
            endpointxdata = startpointxdata+dtbinsize(dt)-1;
            hold on
            if linorlogBGsubtract == 1 %linear
                h=histogram(xdata(startpointxdata:endpointxdata),linspace(0,max(xdata),100),'Normalization','probability');
                histFormat(h);
                Dplot = D1+loc_unc^2/frame_time_dt;
                offsetplot = offset{dt};
                xdatapl = linspace(0,max(xdata),100);
                ydata = (xdatapl.*exp(-(xdatapl.^2)./(4*Dplot*frame_time_dt)))+offsetplot;
                plot(xdatapl,(ydata./sum(ydata).*(size(xdatapl,2)/100)),'k-','LineWidth',2,'DisplayName',num2str(D1))
            elseif linorlogBGsubtract == 2 %Logarithmic
                h=histogram(xdata(startpointxdata:endpointxdata),logspace(log10(minlogpoint),log10(max(xdata)),100),'Normalization','probability');
                histFormat(h);
                set(gca,'XScale','log')
                Dplot = D1+loc_unc^2/frame_time_dt;
                offsetplot = offset{dt};
                xdatapl = logspace(log10(minlogpoint),log10(max(xdata)),100);
                ydatal = (xdatapl.*exp(-(xdatapl.^2)./(4*Dplot*frame_time_dt))).*xdatapl+offsetplot;
                plot(xdatapl,(ydatal./sum(ydatal).*(size(xdatapl,2)/100)),'k-','LineWidth',2,'DisplayName',num2str(D1))
            end
            xlabel('JD (m)')
            ylabel('Occurance')
            legend()
            title(['DEBUG fitting in progress - dt:' num2str(dt)])
        end
        drawnow
    elseif populations == 2
        figure(92);clf(92);
        for dt = 1:size(dtbinsize,2)
            frame_time_dt = frame_time*dt;
            subplot(1,size(dtbinsize,2),dt)
            if dt > 1
                startpointxdata = sum(dtbinsize(1:dt-1))+1;
            else
                startpointxdata = 1;
            end
            endpointxdata = startpointxdata+dtbinsize(dt)-1;
            hold on
            if linorlogBGsubtract == 1 %linear
                h=histogram(xdata(startpointxdata:endpointxdata),linspace(0,max(xdata),100),'Normalization','probability');
                histFormat(h);
                Dplot1 = D1+loc_unc^2/frame_time_dt;
                offset{dt} = (D1+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
                offsetplot = offset{dt}*(1/(1+f_pop_1));
                xdatapl = linspace(0,max(xdata),100);
                ydata1 = (xdatapl.*exp(-(xdatapl.^2)./(4*Dplot1*frame_time_dt)))+offsetplot;
                plot(xdatapl,(ydata1./sum(ydata1).*(size(xdatapl,2)/100)).*(1/(1+f_pop_1)),'b-','LineWidth',2,'DisplayName',num2str(D1))
                
                offset{dt} = (D2+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
                Dplot2 = D2+loc_unc^2/frame_time_dt;
                offsetplot = offset{dt}*(f_pop_1/(1+f_pop_1));
                xdatapl = linspace(0,max(xdata),100);
                ydata2 = (xdatapl.*exp(-(xdatapl.^2)./(4*Dplot2*frame_time_dt)))+offsetplot;
                plot(xdatapl,(ydata2./sum(ydata2).*(size(xdatapl,2)/100)).*(f_pop_1/(1+f_pop_1)),'b-','LineWidth',2,'DisplayName',[num2str(D2) '-' num2str(f_pop_1)])
                
                ydsum = (ydata1./sum(ydata1).*(size(xdatapl,2)/100)).*(1/(1+f_pop_1))+...
                    (ydata2./sum(ydata2).*(size(xdatapl,2)/100)).*(f_pop_1/(1+f_pop_1));
                plot(xdatapl,ydsum,'k-','LineWidth',2);
            elseif linorlogBGsubtract == 2 %Logarithmic
                %                 h=histogram(xdata(startpointxdata:endpointxdata),logspace(log10(minlogpoint),log10(max(xdata)),100),'Normalization','probability');
                %                 histFormat(h);
                %                 set(gca,'XScale','log')
                %                 Dplot = D1+loc_unc^2/frame_time_dt;
                %                 offsetplot = offset{dt};
                %                 xdatapl = logspace(log10(minlogpoint),log10(max(xdata)),100);
                %                 ydatal = (xdatapl.*exp(-(xdatapl.^2)./(4*Dplot*frame_time_dt))).*xdatapl+offsetplot;
                %                 plot(xdatapl,(ydatal./sum(ydatal).*(size(xdatapl,2)/100)),'k-','LineWidth',2,'DisplayName',num2str(D1))
            end
            xlabel('JD (m)')
            ylabel('Occurance')
            legend()
            title(['DEBUG fitting in progress - dt:' num2str(dt)])
        end
        drawnow
    end
end
end