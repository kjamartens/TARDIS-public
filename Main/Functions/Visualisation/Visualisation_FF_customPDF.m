%% Also accepts full appinfo as input!

function fig = Visualisation_FF_customPDF(UIinfo,figID,visInfoCell)
%Differentiate between giving appinfo or visInfoCell
if ~isfield(visInfoCell,'size_dt')
    appinfo = visInfoCell;
    visInfoCell = visInfoCell.Visualisation_FF_outputCell;
end

size_dt = visInfoCell.size_dt;
JDarrSignalCell = visInfoCell.JDarrSignalCell;
xdataplot = visInfoCell.xdataplot;
BGcurve_interp = visInfoCell.BGcurve_interp;
FFparameters = visInfoCell.FFparameters;
fitwithBleach = visInfoCell.fitWithBleach;
populations = visInfoCell.populations;
linorlog = visInfoCell.linorlog;
extraoutput = visInfoCell.extraoutput;
linorlog = 'lin'; %for now? forever?
customPDF = visInfoCell.customPDF;


%% Check if JDarrSignalCell is stored, if not, re-calculate it
if isempty(JDarrSignalCell)
    for dt = appinfo.settingsURDA.dt_arr %Loop over all dt bins wanted
    JDarrSignalCell{dt} = JD_Rel(visInfoCell.poslist,dt,appinfo.settingsURDA.maxdist,0); %Get the relative JDs for every dt bin
    size_dt_new(dt) = size(JDarrSignalCell{dt},1); %Get the number of JD values for every dt bin
    
    %Cleanup of zeros
    JDarrSignalCell{dt}(JDarrSignalCell{dt}==0) = [];

    %Check if OK with size_dt
    if size(size_dt_new) ~= size(size_dt)
        disp('OHNO')
    end
    end
end

%%
if fitwithBleach
    arr_parameters_bleach = FFparameters; %Transform the parameters from cell to array
        indexFirstBG = size(FFparameters,2)-1;
        indexFitParams = [1:(size(arr_parameters_bleach,2)-2)];
        BGfDt(1) = arr_parameters_bleach(indexFirstBG);
        bleachtimeframes = arr_parameters_bleach(end)./extraoutput.frame_time;
        %Now we get the decay based on the bleach time only - this still
        %needs to be corrected for nr of JDs and BG fraction at dt = 1
        BGfDt_decay = multiExpFitTARDISBleach(bleachtimeframes,[1:size(size_dt,2)],[1:size(size_dt,2)]);
        %First corrected for nr of JDs
        for i = 1:size(JDarrSignalCell,2); nrlinks(i) = size(JDarrSignalCell{1,i},1); end
        BGfDt_decay = BGfDt_decay./nrlinks';
        BGfDt_decay = BGfDt_decay./max(BGfDt_decay);
        %Now using this to calculate the BGfractions
        SignalFractions = (1-BGfDt(1)).*BGfDt_decay;
        BGfDt = 1-SignalFractions;
        %Restore varargin
        FFparameters = [arr_parameters_bleach(indexFitParams)'; fliplr(BGfDt')'];
else
    indexFirstBG = size(FFparameters,2)-size(size_dt,2)+1;
end
%% Prepare plot area
maxheight = 0;
if isempty(UIinfo)
    fig = figure(5);clf(5);
    plot([0,0]);
elseif UIinfo == 'hidden'
    fig = figure('visible','off'); clf();
else
    eval(['UIinfo.' figID '.AutoResizeChildren = ''off'';']);
end
% Make the plots
maxdtbins = max(size(size_dt));
maxxpos = [];
for dt = 1:maxdtbins %Make subfigures for all dt
    if isempty(figID)%Make the subplot
        ax = subplot(5,maxdtbins,[dt, dt+maxdtbins, dt+2*maxdtbins]);
    else
        eval(['ax = subplot(5,maxdtbins,[dt, dt+maxdtbins, dt+2*maxdtbins],''Parent'',UIinfo.' figID ');']);
    end

    startpos = sum(size_dt(1:dt))-size_dt(1)+1; %Get the startpos of the arrays based on dt size
    endpos = startpos+size_dt(dt)-1; %Get the endpos of the arrays based on dt size
    if linorlog == 'lin'
        %Plot histogram of data
        xdataplot2 = linspace(min(xdataplot(startpos:endpos)),max(xdataplot(startpos:endpos)),size(xdataplot(startpos:endpos),1)+0);
        h = histogram(ax,JDarrSignalCell{dt},xdataplot2,'Normalization','probability','DisplayName','Data');
        histArea = sum(h.Values); %Find the histogram area
        histData = h.Values;
        histFormat(h);
        hold(ax,'on');

        %Plot the BG curve
        BGydata = BGcurve_interp(:,2)./sum(BGcurve_interp(:,2)); %Plot noise data - normalized
        BGydata = BGydata.*FFparameters(indexFirstBG+maxdtbins-(dt-1)-1); %Correct for BG percentage found
        plot(ax,BGcurve_interp(:,1),BGydata.*histArea,'r-','LineWidth',1,'DisplayName','Background');

        %Plot the signal curve
        signalxdata = xdataplot(startpos:endpos);
        signalydata = extraoutput.main{dt}{2};
        signalydata = signalydata./sum(signalydata);
        signalydata = signalydata.*(1-FFparameters(indexFirstBG+maxdtbins-(dt-1)-1)); %Correct for signal percentage found
        plot(ax,signalxdata,signalydata,'b-','LineWidth',1,'DisplayName','Signal');

        %Plot summed populations and BG
        sumdata = BGydata.*histArea+signalydata; 
        plot(ax,xdataplot(startpos:endpos),sumdata,'k-','LineWidth',1,'DisplayName','Overall fit'); %Draw it

        %Calculate residuals
        %Fitted sumdata is 101 values at 0, 5, 10, etc.
        %Bins are 100 bins with midpoints 2.5, 7.5, etc.
        % Solution: Plot everything at midpoints of bin, interpolate
        % fitted sumdata between the points.
        
        %Get interpolated points of sumdata
        interpsumdata = interp(sumdata,2);
        interpsumdata = interpsumdata(2:2:end-2);
        %Get interpolated points of xdataplot
        interpxdataplot = interp(xdataplot(startpos:endpos),2);
        interpxdataplot = interpxdataplot(2:2:end-2);

        residuals{dt}(:,1) = interpxdataplot;
        residuals{dt}(:,2) = h.Values'-interpsumdata;
%     elseif linorlog == 'log'
%         %% Single-dual pop, log
%         logxdata = extraoutput.logvis.xdataplotlog(startpos:endpos);
%         %Plot histogram of data
%         %Add an extra position in logxdata (i.e. go from 100 datapoints to
%         %101 datapoints - to get 100 bins with histogram())
%         logxdata2 = logspace(log10(min(logxdata)),log10(max(logxdata)),size(logxdata,1)+1)';
%         h = histogram(ax,JDarrSignalCell{dt},logxdata2,'Normalization','probability','DisplayName','Data');
%         histArea = sum(h.Values); %Find the histogram area
%         histData = h.Values;
%         histFormat(h);
%         %Get the middle of the logbins (taking half a bin into account)
%         %These should be the x-positions of the plotted curves later.
%         histMidLogBins = 10.^(log10(h.BinEdges(1:end-1))+(log10(h.BinEdges(2))-log10(h.BinEdges(1)))/2);
% 
%         hold(ax,'on');
%         set(ax,'XScale','log')
%         %Plot the BG curve
%         BGydata = (extraoutput.logvis.BGcurve_interplog(:,2)./sum(extraoutput.logvis.BGcurve_interplog(:,2))); %Plot noise data - normalized
%         BGydata = BGydata(1:end).*FFparameters(end-(dt-1)); %Correct for BG percentage found
%         %X positions to plot data on should be the mid positions of the
%         %logbins
%         logxdata = histMidLogBins;
%         plot(ax,logxdata,BGydata.*histArea,'r-','LineWidth',1,'DisplayName','Background');
% 
%         %Plot the populations
%         for p = 1:populations
%             pxdata = extraoutput.logvis.xdataplotlog(startpos:endpos); %Get xdata
% 
%             %Generate ydata
%             %pydata = pxdata.*pxdata.*exp((-pxdata.^2)/(4*(FFparameters(p)*1e-12+extraoutput.loc_unc.^2./extraoutput.frame_time)*(extraoutput.frame_time*dt)));
%             %Correcting for strobo illumination here                
%             FFparameters(p) = FFparameters(p)+FFparameters(p)*((extraoutput.strobo_frame_time*(1/6))/(extraoutput.frame_time*4*dt));
% 
%             pydata = pxdata.*pxdata.*exp((-pxdata.^2)/(4*(FFparameters(p)*1e-12+extraoutput.loc_unc.^2./(extraoutput.frame_time))*(extraoutput.frame_time*dt)));
%             pydata = pydata./sum(pydata); %Normalize ydata (same normalization as for BG curve)
%             %Correct for population ratio and BG percentage found
%             startparam_posratio = populations+1; %Get the startposition in the parameters list for lookup
%             endparam_posratio = startparam_posratio+populations-2; %Get the endposition in the parameters list for lookup
%             %Exception for 1 population, since it's always a 'non-existant'
%             %ratio between populations
%             if populations == 1
%                 ratio_pop = 1;
%             else % in the case of having more than 1 population
%                 paramval_posratio = startparam_posratio+(p-2); %Get the value in the position array where the posratio is calculated
%                 if p == 1 %if it's the ratio for the first population
%                     ratio_pop = (1/(1+sum(FFparameters(startparam_posratio:endparam_posratio))));
%                 else %else if it's the ratio for the second popultion - NOTE that it should be implemented for a third population
%                     ratio_pop = (FFparameters(paramval_posratio)/(1+sum(FFparameters(startparam_posratio:endparam_posratio))));
%                 end
%             end
%             %Correct for population and BG ratio - BG ratio is dt-dependant
%             pydata = pydata.*ratio_pop.*(1-FFparameters(end-(dt-1)));
%             plot(ax,logxdata,pydata.*histArea,'b-','LineWidth',1,'DisplayName',['Pop. ' num2str(p)]); %Plot data and change layout
%             pydata_cell{p} = pydata.*histArea; %Save the population data in a cell array for later drawing the summed population
%             if dt == 1
% %                     xlimmaxv = logxdata(pydata==max(pydata));
%                 ylimv = max(pydata_cell{p})*1.8;
%                 xlimv = 2e-6;
% %                     disp(['Suggested Y lim: ' num2str(max(pydata_cell{p})*1.33)]);
%             end
%         end
%         %Plot summed populations and BG
%         sumdata = BGydata.*histArea; %Start creating summed data with BG-ydata, normalized to histArea
%         for p = 1:populations
%             try
%                 sumdata = sumdata+pydata_cell{p}; %Add populations-ydata. Note that these are already normalized to histArea.
%             catch
%                 keyboard
%             end
%         end
%         plot(ax,logxdata,sumdata,'k-','LineWidth',1,'DisplayName','Overall fit'); %Draw it
%         set(ax,'XScale','log')
% 
%         %Calculate residuals
%         residuals{dt}(:,1) = logxdata;
%         residuals{dt}(:,2) = histData'-sumdata;
    end
    xlabel(ax,'JD (m)')

    if dt == 1
        ylabel(ax,'Occurance')
    else
        set(ax,'YTick',[]);
    end
    title(ax,['TARDIS fit at \Deltat = ' num2str(dt)]);
    if dt == maxdtbins
        %         lgd = legend(ax);
        %         lgd.FontSize = 9;
        %         title(lgd,'URDA Legend')
        %         lgd.Position(1) = lgd.Position(1)-0.575;
    end

    maxheightnew = max(max(sumdata),max(histData))*1.1;
    if maxheightnew > maxheight
        maxheight = maxheightnew;
    end
%     maxheight = ylimv;
end
%% Alter the axis to have same max height
for dt = 1:maxdtbins %Make subfigures for all dt
    if isempty(figID)%Make the subplot
        ax = subplot(5,maxdtbins,[dt, dt+maxdtbins, dt+2*maxdtbins]);
    else
        eval(['ax = subplot(5,maxdtbins,[dt, dt+maxdtbins, dt+2*maxdtbins],''Parent'',UIinfo.' figID ');']);
    end
    if linorlog == 'lin'
        axis(ax,[0 BGcurve_interp(end,1) 0 maxheight]);
    else
        axis(ax,[extraoutput.logvis.xdataplotlog(1) extraoutput.logvis.xdataplotlog(end) 0 maxheight]);
    end
end
%% Create residuals plots
for dt = 1:maxdtbins %Make subfigures for all dt
    if isempty(figID)%Make the subplot
        ax = subplot(5,maxdtbins,[dt+4*maxdtbins]);
    else
        eval(['ax = subplot(5,maxdtbins,[dt+4*maxdtbins],''Parent'',UIinfo.' figID ');']);
    end
    %Change pos/height a little
    ax.Position(2) = ax.Position(2) + ax.Position(4)*0.2; %Shift up
    ax.Position(4) = ax.Position(4)*1.4; %Increase height

    plot(ax,residuals{dt}(1:end-1,1),residuals{dt}(1:end-1,2),'k-')
    minmaxy = 1.1*max([abs(max(residuals{dt}(1:end-1,2))), abs(min(residuals{dt}(1:end-1,2)))]);

    xlabel(ax,'JD (m)')
    if linorlog == 'lin'
        axis(ax,[0 BGcurve_interp(end,1) -minmaxy minmaxy]);
    else
        axis(ax,[extraoutput.logvis.xdataplotlog(1) extraoutput.logvis.xdataplotlog(end) -minmaxy minmaxy]);
        set(ax,'XScale','log');
    end
   
    %Ylabel at dt=1
    if dt == 1
        ylabel(ax,'Residual')
    end
end
drawnow; 
end