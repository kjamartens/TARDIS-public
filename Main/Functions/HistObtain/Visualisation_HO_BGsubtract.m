function fig = Visualisation_HO_BGsubtract(UIinfo,figID,bgbinningnr,dt_arr,JDarrBG,maxdist,minlogpoint,linorlogBGsubtract,bgarr,JDarrSignalCell,...
    signalCurve_interp,stepsizearraytrue,bgratios,startpointBG,app)
% keyboard
if ~isempty(app)
    nrpopulations = app.Data.settingsURDA.populations;
    frame_time = app.Data.settingsURDA.frame_time;
end
%% Visualisation outside MLE_BG_subtraction_HistObtain_function
maxhheight = 0;
nrbinsvis = bgbinningnr+1;
%Create figure
if isempty(UIinfo)
fig = figure(1);clf(1);
plot([0,0]); 
gca = fig.CurrentAxes;
elseif UIinfo == 'hidden'
    fig = figure('visible','off'); clf();
else
    eval(['UIinfo.' figID '.AutoResizeChildren = ''off'';']);
end
%Make subplots based on how many dt-bins I have
%% BG-only plot:
if isempty(figID)
    ax = subplot(2,size(dt_arr,2)+1,1);
else
    eval(['ax = subplot(2,size(dt_arr,2)+1,1,''Parent'',UIinfo.' figID ');']);
end
if linorlogBGsubtract == 'lin'
    h=histogram(ax,JDarrBG,linspace(0,maxdist,nrbinsvis),'Normalization','probability');
elseif linorlogBGsubtract == 'log'
    h=histogram(ax,JDarrBG,logspace(log10(minlogpoint),log10(maxdist),nrbinsvis),'Normalization','probability');
end
histFormat(h);
hold(ax,'on')
if linorlogBGsubtract == 'lin'
    plotydata = [0; bgarr(:,2).*((bgbinningnr+1)/nrbinsvis)];
    plot(ax,[0; bgarr(:,1)],plotydata,'k-');
elseif linorlogBGsubtract == 'log'
    plotydata = ([bgarr(:,2).*((bgbinningnr+1)/nrbinsvis)]);
    plot(ax,[bgarr(:,1)],plotydata,'k-');
end
if max(plotydata)>maxhheight; maxhheight = max(plotydata);end
plot(ax,[startpointBG,startpointBG],[0,maxhheight*90],'k--')
title(ax,'Background info')

%% BG+signal plots on dt
for dt = dt_arr
    
    if isempty(figID)
        ax = subplot(2,size(dt_arr,2)+1,1+dt);
    else
    eval(['ax = subplot(2,size(dt_arr,2)+1,1+dt,''Parent'',UIinfo.' figID ');']);
    end

    if linorlogBGsubtract == 'lin'
        h=histogram(ax,JDarrSignalCell{dt},linspace(0,maxdist,nrbinsvis),'Normalization','probability');
    elseif linorlogBGsubtract == 'log'
        if nrpopulations > 0
            h=histogram(ax,JDarrSignalCell{dt},logspace(log10(minlogpoint),log10(maxdist),nrbinsvis),'Normalization','probability');
        elseif nrpopulations == 0 %anaDDA
            plotx = (JDarrSignalCell{dt}.*1e6).^2./(4*frame_time*dt);
            minlogpointD = (minlogpoint.*1e6).^2./(4*frame_time*dt);
            maxdistD = (maxdist.*1e6).^2./(4*frame_time*dt);
            h=histogram(ax,plotx,logspace(log10(minlogpointD),log10(maxdistD),nrbinsvis),'Normalization','probability');
        end
    end
    histFormat(h);
    hold(ax,'on')
    if linorlogBGsubtract == 'lin'
        plot(ax,[0; bgarr(:,1)],[0; bgarr(:,2)*bgratios(dt)],'r-');
        plotydata = [0; signalCurve_interp{dt}(:,2).*((bgbinningnr+1)/nrbinsvis)];
        plot(ax,[0; signalCurve_interp{dt}(:,1)],plotydata,'k-');
    elseif linorlogBGsubtract == 'log'
        if nrpopulations > 0
            plot(ax,[bgarr(:,1)],[bgarr(:,2)*bgratios(dt)],'r-');
            plotydata = [signalCurve_interp{dt}(:,2).*((bgbinningnr+1)/nrbinsvis)];
            plot(ax,[signalCurve_interp{dt}(:,1)],plotydata,'k-');
        elseif nrpopulations == 0 %anaDDA
            bgarr1D = (bgarr(:,1).*1e6).^2./(4*frame_time*dt);
            signalcurveinterp1 = (signalCurve_interp{dt}(:,1).*1e6).^2./(4*frame_time*dt);
            plot(ax,[bgarr1D],[bgarr(:,2)*bgratios(dt)],'r-');
            plotydata = [signalCurve_interp{dt}(:,2).*((bgbinningnr+1)/nrbinsvis)];
            plot(ax,[signalcurveinterp1],plotydata,'k-');
            startpointBG_D = (startpointBG.*1e6).^2./(4*frame_time*dt);
            plot(ax,[startpointBG_D,startpointBG_D],[0,maxhheight*90],'k--')
        end
    end
    plot(ax,[startpointBG,startpointBG],[0,maxhheight*90],'k--')
    if max(plotydata)>maxhheight; maxhheight = max(plotydata);end
    title(ax,['Info of \Deltat=' num2str(dt)])
end

%% Set axis of top row
for i = 1:size(dt_arr,2)+1
    
    if isempty(figID)
        ax = subplot(2,size(dt_arr,2)+1,i);
    else
    eval(['ax = subplot(2,size(dt_arr,2)+1,i,''Parent'',UIinfo.' figID ');']);
    end
    
    if linorlogBGsubtract == 'lin'
        axis(ax,[0 maxdist 0 maxhheight*1.1])
    xlabel(ax,'JD (m)')
    elseif linorlogBGsubtract == 'log'
        if nrpopulations > 0 || i == 1
            axis(ax,[minlogpoint maxdist 0 maxhheight*1.1]);
            set(ax,'XScale','log')
            xlabel(ax,'JD (m)')
        elseif nrpopulations == 0
            maxdistD = (maxdist.*1e6).^2./(4*frame_time*(i-1));
            axis(ax,[minlogpointD maxdistD 0 maxhheight*1.1])
            set(ax,'XScale','log')
            xlabel(ax,'D* (\mum^2/s)')
        end
    end
    if i > 1
%         set(gca,'YTick',[]);
        set(ax,'YTick',[]);
    else
        ylabel(ax,'Occurance');
    end
end

%% Signal only plots
multfactor = 1;
for dt = dt_arr
    if isempty(figID)
        ax = subplot(2,size(dt_arr,2)+1,size(dt_arr,2)+1+dt+1);
    else
    eval(['ax = subplot(2,size(dt_arr,2)+1,size(dt_arr,2)+1+dt+1,''Parent'',UIinfo.' figID ');']);
    end
    
    if linorlogBGsubtract == 'lin'
        h=histogram(ax,stepsizearraytrue{dt},linspace(0,maxdist,nrbinsvis*multfactor),'Normalization','probability');
    elseif linorlogBGsubtract == 'log'
        if nrpopulations > 0
            h=histogram(ax,stepsizearraytrue{dt},logspace(log10(minlogpoint),log10(maxdist),nrbinsvis),'Normalization','probability');
        elseif nrpopulations == 0
            stepsizearraytruedataD = (stepsizearraytrue{dt}.*1e6).^2./(4*frame_time*dt);
            startpointBG_D = (startpointBG.*1e6).^2./(4*frame_time*dt);
            plot(ax,[startpointBG_D,startpointBG_D],[0,maxhheight*90],'k--')
            hold(ax,'on');
            maxdistD = (maxdist.*1e6).^2./(4*frame_time*dt);
            h=histogram(ax,stepsizearraytruedataD,logspace(log10(minlogpointD),log10(maxdistD),nrbinsvis),'Normalization','probability');
        end
    end
    h.FaceColor = [0.5 0.5 0.5];
    h.EdgeColor = 'none';
    hold(ax,'on');
    ratio = bgratios(dt);
    if linorlogBGsubtract == 'lin'
        datax = [0; bgarr(:,1)];
        datay = [0; signalCurve_interp{dt}(:,2) - bgarr(:,2)*ratio];
    elseif linorlogBGsubtract == 'log'
        if nrpopulations > 0
            datax = [bgarr(:,1)];
        elseif nrpopulations == 0
            datax = (bgarr(:,1).*1e6).^2./(4*frame_time*dt);
        end
        datay = [signalCurve_interp{dt}(:,2) - bgarr(:,2)*ratio];
    end
    plotdatay{dt} = datay./(1-ratio).*((bgbinningnr+1)/(nrbinsvis*multfactor));
    plotdatay{dt} = plotdatay{dt}./sum(plotdatay{dt});
    plot(ax,datax,plotdatay{dt},'k-')
    plot(ax,datax,zeros(size(datax)),'k:')
    plot(ax,[startpointBG,startpointBG],[-100,max(max(cell2mat(plotdatay)))*2],'k--')
    title(ax,['Data minus BG, \Deltat=' num2str(dt)]);
end

%% HO fit show
try
if app.Data.settingsURDA.performestimationfit
    if nrpopulations == 0
        size_dt = app.Data.anaDDAvisInfoHO.size_dt_HO;
        for dt = 1:size(app.Data.settingsURDA.dt_arr,2)
            if isempty(figID)%Make the subplot
                ax = subplot(2,size(app.Data.settingsURDA.dt_arr,2)+1,size(app.Data.settingsURDA.dt_arr,2)+1+dt+1);
            else
                eval(['ax = subplot(2,size(app.Data.settingsURDA.dt_arr,2)+1,size(app.Data.settingsURDA.dt_arr,2)+1+dt+1,''Parent'',UIinfo.' figID ');']);
            end
%             ax = subplot(2,size(app.Data.settingsURDA.dt_arr,2)+1,size(app.Data.settingsURDA.dt_arr,2)+1+dt+1,'Parent',app.RightTopPanel);
            hold(ax,'on');
            
            set(ax,'XScale','log')
            startpos = sum(size_dt(1:dt))-size_dt(dt)+1; %Get the startpos of the arrays based on dt size
            endpos = startpos+size_dt(dt)-1; %Get the endpos of the arrays based on dt size
            
            %The only input we have is from some archaic visualisation which
            %gives a bunch of simulated xvalues and corresponding yvalues.
            %They're not neatly organised
            xdatarandvalues = app.Data.anaDDAvisInfoHO.DlistHO(startpos:endpos);
            ydatarandvalues = app.Data.anaDDAvisInfoHO.visual{1,dt}{1,2}{1,dt};
            
            %Sort these values
            [~,srt] = sort(xdatarandvalues);
            xdatarandvalues = xdatarandvalues(srt);
            ydatarandvalues = ydatarandvalues(srt);
            
            %Now get the xvalues of D-ified-bgarr(:,1) as interpolated values of ydata
            bgarrxD = (bgarr(:,1).*1e6).^2./(4*app.Data.settingsURDA.frame_time*dt);
            t = interp1(xdatarandvalues,ydatarandvalues,bgarrxD);
            %Cleanup and normalize
            t(isnan(t)) = 0;
            t = t./sum(t);
            %Plot
            plot(ax,bgarrxD,t,'b-')
        end
        
    else %1 or 2 populations - plot the estimated population
        for dt = dt_arr
        if isempty(figID)%Make the subplot
            ax = subplot(2,size(app.Data.settingsURDA.dt_arr,2)+1,size(app.Data.settingsURDA.dt_arr,2)+1+dt+1);
        else
            eval(['ax = subplot(2,size(app.Data.settingsURDA.dt_arr,2)+1,size(app.Data.settingsURDA.dt_arr,2)+1+dt+1,''Parent'',UIinfo.' figID ');']);
        end
%         ax = subplot(2,size(app.Data.settingsURDA.dt_arr,2)+1,size(app.Data.settingsURDA.dt_arr,2)+1+dt+1,'Parent',app.RightTopPanel);
        hold(ax,'on');
        if nrpopulations == 1
            ratiov(1) = 1;
        elseif nrpopulations == 2
            ratiov(2) = app.Data.paramEsts(3)/(1+app.Data.paramEsts(3));
            ratiov(1) = 1-ratiov(2);
        end
        for p = 1:nrpopulations
            %Get main D
            DHO{p} = app.Data.paramEsts(p)*1e-12;
            sigma = app.Data.settingsURDA.loc_unc;
            dtval = app.Data.settingsURDA.frame_time*dt;
            if strcmp(app.Data.settingsURDA.linorlogBGsubtract,'lin')
                xdata = linspace(0,app.Data.settingsURDA.maxdist,app.Data.settingsURDA.bgbinningnr);
                ydata = xdata.*exp((-xdata.^2)./((4*DHO{p}+sigma.^2/dtval)*dtval));
                %Normalise and include ratio
                ydata = ydata./sum(ydata).*ratiov(p);
            else
                xdata = logspace(log10(minlogpoint),log10(app.Data.settingsURDA.maxdist),app.Data.settingsURDA.bgbinningnr);
                ydata = xdata.*(xdata.*exp((-xdata.^2)./((4*DHO{p}+sigma.^2/dtval)*dtval)));
                %Normalise and include ratio
                ydata = ydata./sum(ydata).*ratiov(p);
            end
            plot(ax,xdata,ydata,'b-.')
        end
        end
    end
end
%% Set axis nice and same, so loop back over the same subplots as above
for dt = dt_arr
    if isempty(figID)
        ax = subplot(2,size(dt_arr,2)+1,size(dt_arr,2)+1+dt+1);
    else
    eval(['ax = subplot(2,size(dt_arr,2)+1,size(dt_arr,2)+1+dt+1,''Parent'',UIinfo.' figID ');']);
    end
    if linorlogBGsubtract == 'lin'
        axis(ax,[0 maxdist min(min(cell2mat(plotdatay)))*1.1 max(max(cell2mat(plotdatay)))*1.1])
        xlabel(ax,'JD (m)')
    elseif linorlogBGsubtract == 'log'
        if nrpopulations > 0
            %Clause for when for whatever reason the limits are wrong -
            %happened once...
            yminlim = min(min(cell2mat(plotdatay)))*1.1;
            ymaxlim = max(max(cell2mat(plotdatay)))*1.1;
            if isnan(yminlim)
                yminlim = 0;
            end
            if isnan(ymaxlim)
                ymaxlim = inf;
            end
            axis(ax,[minlogpoint maxdist yminlim ymaxlim])
            xlabel(ax,'JD (m)')
        elseif nrpopulations == 0
            maxdistD = (maxdist.*1e6).^2./(4*app.Data.settingsURDA.frame_time*dt);
            axis(ax,[minlogpointD maxdistD min(min(cell2mat(plotdatay)))*1.1 max(max(cell2mat(plotdatay)))*1.1])
            xlabel(ax,'D* (\mum^2/s)')
        end
        set(ax,'XScale','log')
    end
    if dt > 1
        set(ax,'YTick',[]);
    else
        ylabel(ax,'Occurance');
    end
end
catch
    disp('Skipped a part here - probably OK if running from code')
end
drawnow
end
