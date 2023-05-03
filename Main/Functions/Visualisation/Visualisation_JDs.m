function fig = Visualisation_JDs(UIinfo,figID,JDcelldata,Datacell)
%% Prepare plot area
if isempty(UIinfo)
    fig = figure(7);clf(7);
    plot([0,0]);
elseif UIinfo == 'hidden'
    fig = figure('visible','off'); clf();
else
    eval(['UIinfo.' figID '.AutoResizeChildren = ''off'';']);
end

%% Make the plots
maxdtbins = size(JDcelldata,2);
nrrows = 2; %once linear, once log
nrbinshist = 50;
%% Make the linear plots
maxheight = 0;
for dt = 1:maxdtbins %Make subfigures for all dt    
    if isempty(figID)%Make the subplot
        ax = subplot(2,maxdtbins,dt);
    else
        eval(['ax = subplot(2,maxdtbins,dt,''Parent'',UIinfo.' figID ');']);
    end
    h = histogram(ax,JDcelldata{dt},linspace(0,Datacell.Data.settingsURDA.maxdist,nrbinshist),'Normalization','probability');
    histFormat(h);
    newmaxheight = max(h.Values);
    if newmaxheight > maxheight
        maxheight = newmaxheight;
    end
end
%Set ax limits
for dt = 1:maxdtbins
    if isempty(figID)%Make the subplot
        ax = subplot(2,maxdtbins,dt);
    else
        eval(['ax = subplot(2,maxdtbins,dt,''Parent'',UIinfo.' figID ');']);
    end
    axis(ax,[0 Datacell.Data.settingsURDA.maxdist 0 maxheight*1.1]);
    title(ax,['Approx. JDs at \Deltat = ' num2str(dt)]);
    xlabel(ax,'JD (m)')
    if dt == 1
        ylabel(ax,'Probability')
    end
end
%% Make the log plots
maxheight = 0;
for dt = 1:maxdtbins %Make subfigures for all dt    
    if isempty(figID)%Make the subplot
        ax = subplot(2,maxdtbins,dt+maxdtbins);
    else
        eval(['ax = subplot(2,maxdtbins,dt+maxdtbins,''Parent'',UIinfo.' figID ');']);
    end
    h = histogram(ax,JDcelldata{dt},logspace(log10(Datacell.Data.settingsURDA.minlogpoint),log10(Datacell.Data.settingsURDA.maxdist),nrbinshist),'Normalization','probability');
    histFormat(h);
    set(ax,'XScale','log')
    newmaxheight = max(h.Values);
    if newmaxheight > maxheight
        maxheight = newmaxheight;
    end
end
for dt = 1:maxdtbins
    if isempty(figID)%Make the subplot
        ax = subplot(2,maxdtbins,dt+maxdtbins);
    else
        eval(['ax = subplot(2,maxdtbins,dt+maxdtbins,''Parent'',UIinfo.' figID ');']);
    end
    axis(ax,[Datacell.Data.settingsURDA.minlogpoint Datacell.Data.settingsURDA.maxdist 0 maxheight*1.1]);
    title(ax,['Approx. JDs at \Deltat = ' num2str(dt)]);
    xlabel(ax,'JD (m)')
    if dt == 1
        ylabel(ax,'Probability')
    end
end
%%
end