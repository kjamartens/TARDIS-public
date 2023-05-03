%% Create 'reconstituted' JD array
% Create a series of arrays of JD values that correspond to basically the 
% 'BG-subtracted' information in TARDIS. The basic idea is to remove the BG
% data from the dt bins and then sub-/over-sample the resulting JDs in a
% reconstituted array
%---------------------------------------------------------
% Required inputs
% JDarrSignalCell:      Cell containing the JD info for every DT
% bgratios:             Ratios of BG that is found (0-1) for every dt
% linorlogBGsubtract:   'lin' or 'log', depending on a linear or
%                       logarithmic BG subtract method
% histInfo:             cell for every dt that contains [x,y]-data of the
%                       signal curve - created via MLE_BG_subtraction_HistObtain_function
% debug:                Boolean to get some extra debug info
% maxdist:              Max x distance for the plot (m)
% minlogpoint:          Minimum x value in case log value is used
% dt_arr:               Array of used dt values. Always increasing from 1
%                       upwards
% bgbinningnr:          Number of bins used for the BG correction
%
% Output
% reconstit_arr:        Cells of the JD reconstituted array for every dt
% JDonlydata:           Placing the reconstit_arr below each other, basically
% truthoffsetpartial:   Absolute value of the baseline that is added for
%                       every dt
% baseline_to_be_added: y-Offset of the whole JD array so that there are no
%                       negative values
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [JDonlydata, reconstit_arr, truthoffsetpartial, baseline_to_be_added] = createReconstitutedJDarray(JDarrSignalCell,bgratios,linorlogBGsubtract,histInfo,debug,maxdist,minlogpoint,dt_arr,bgbinningnr)
%Get number of entries in total - note: 1k rather than 10k at HO
maxJDvaluesWantedForJDarray = max(ceil(cellfun('size',JDarrSignalCell,1).*(1-bgratios)));
totalreconstit_arr_entries = max(50000,ceil(1.1*(1+maxJDvaluesWantedForJDarray))); % used to be 10000 % This is a sort of minimum? It's slightly higher than this because of removal of negatives
%Enable or disable having a 'baseline' for every bin. If this is disabled,
%results are likely biased
baseline_per_bin = true;
endpoint = maxdist;
%Find Edges of the histogram
if linorlogBGsubtract == 'lin'
    startpoint = 0;
    halfwidthsize = (histInfo{1}(2,1)-histInfo{1}(1,1))/2;
    if (halfwidthsize-histInfo{1}(1,1)) > 1e-16
        disp('This is going to fail! Probably change code')
    end
    %For linear, simply interpolate between start and endpoint
    edges = [startpoint:halfwidthsize*2:endpoint];
elseif linorlogBGsubtract == 'log'
    %For log, interpolate in log-space
    loghalfwidthsize = (log10(histInfo{1}(2,1))-log10(histInfo{1}(1,1)))/2;
    edges = [log10(minlogpoint):loghalfwidthsize*2:log10(endpoint)];
    edges = 10.^edges;
end
%Loop over the dt bins
for dt = dt_arr
    %Create previous and next interpoints, which dictate the shape of how
    %the random points in the reconstit_arr should be arranged.
    %F.e. a bin that's part of a quickly rising bin, shouldn't be a
    %'block', but more like a slope.
    %Idea is to get nr of JDs in the reconstit_arr from the histInfo, but the
    %shape from prev and next-interp.
    for i = 1:size(edges,2)-1
        if i == 1
            previnterp(i) = 0;
        elseif i == size(edges,2)-1
            nextinterp(i) = 0;
        end
        
        if i > 1
            previnterp(i) = 0.5*(histInfo{dt}(i,2)+histInfo{dt}(i-1,2));
        end
        if i < size(edges,2)-1
            nextinterp(i) = 0.5*(histInfo{dt}(i,2)+histInfo{dt}(i+1,2));
        end
    end
    %Figure for debugging:
    if debug
        if linorlogBGsubtract=='lin'
            figure(90);clf(90);
            hold on
            plot(histInfo{dt}(:,1),histInfo{dt}(:,2),'k-x');
            plot(histInfo{dt}(:,1)-halfwidthsize,previnterp,'r-o');
            plot(histInfo{dt}(:,1)+halfwidthsize,nextinterp,'r-x');
            title('DEBUG on interpolation of non-static data formation')
        end
    end
    %Make empty reconstit_array
    reconstit_arr{dt} = zeros(totalreconstit_arr_entries,1);
    counter = 1;
    baseline_to_be_added{dt} = 0;
    if baseline_per_bin
        for i = 1:size(edges,2)-1
            %find nr of entries per bin for every bin
            nrentries_bin_check(i) = floor(histInfo{dt}(i,2)*totalreconstit_arr_entries);
        end
        %find the baseline that should be added - should be a positive nr 1
        %higher than the smallest negative bin in the histogram
        baseline_to_be_added{dt} = max(0,min(nrentries_bin_check)*-1+1);
    end
    for i = 1:size(edges,2)-1
        %For every bin, get the number of entries based on the size of the
        %bin and the baseline
        if ~baseline_per_bin
            nrentries_bin = max(round(histInfo{dt}(i,2)*totalreconstit_arr_entries),0);
        else
            nrentries_bin = round(histInfo{dt}(i,2)*totalreconstit_arr_entries)+baseline_to_be_added{dt};
        end
        if linorlogBGsubtract == 'lin'
            %Main idea: split up bin in X sub-bins (f.e. 100), with relative
            %probabilities of being chosen from previnterp and nextinterp.
            %Then, reconstit_arr is filled taking care of these
            %probabilities
            
            %This is called 'FancySubBins', and is set to be used for now.
            %Old method is retained
            %Honestly, testing indicates that there's basically no
            %difference anyhow
            %Thus, we skip it for the time being
            fancysubbins = 0;
            movsumnextinterp = movsum(nextinterp,2);
            if min(movsumnextinterp(1:end-1)) == 0
                fancysubbins = 0;
            end
            if fancysubbins
                %Create the nr of subbins
                subbins = 100;
                %Get a probability of a localization being in the subbins
                subbinprob = linspace(previnterp(i),nextinterp(i),subbins);
                if sum(subbinprob) > 0
                    subbinprob = subbinprob./sum(subbinprob);
                end
                subbinprobcum = cumsum(subbinprob);
                tarr = zeros(nrentries_bin,1);
                randarr = rand(nrentries_bin,1);
                %Loop over all entries, choose a subbin based on the
                %cumulative probabilities, and throw a random value in this
                %subbin
                for cntr = 1:nrentries_bin
                    try
                        chosenbin = find((subbinprobcum-randarr(cntr))>0,1,'first');
                        tarr(cntr) = edges(i)+((halfwidthsize*2)/subbins)*(chosenbin-1)+rand(1,1)*(halfwidthsize*2)/subbins;
                    catch
                        keyboard
                    end
                end
                %Update reconstit_arr
                reconstit_arr{dt}(counter:counter+nrentries_bin-1) = tarr;
            else
                %Old idea: choosing nrentries_bin random values within this bin
                %without further consideration
                reconstit_arr{dt}(counter:counter+nrentries_bin-1) = rand(nrentries_bin,1)*(halfwidthsize*2)+edges(i);
            end
        elseif linorlogBGsubtract == 'log'
            %Currently, for logarithmic BG subtraction, only the 'Simple'
            %method is implemented
            reconstit_arr{dt}(counter:counter+nrentries_bin-1) = 10.^((rand(nrentries_bin,1)*(loghalfwidthsize*2)+log10(edges(i))));
        end
        counter = counter + nrentries_bin;
    end
    %remove zero-entries
    reconstit_arr{dt}(reconstit_arr{dt}==0) = [];
    % Visualisation of blockyness due to binsize for debugging
    if debug
        if linorlogBGsubtract=='lin'
            figure(91);clf(91);
            mf = 4;
            ratio_baseline_to_data{dt} = (baseline_to_be_added{dt}*(size(edges,2)-1))/size(reconstit_arr{dt},1);
            histogram(reconstit_arr{1},linspace(0,maxdist,mf*(bgbinningnr+1)),'Normalization','probability')
            hold on
            yd = histInfo{1}(:,2)*(1-ratio_baseline_to_data{1})+ones(size(histInfo{1}(:,2)))./bgbinningnr*ratio_baseline_to_data{1};
            plot(histInfo{1}(:,1),(yd./sum(yd))/mf,'r-x');
            title('DEBUG on visualization of blockyness due to binsize')
        end
    end
    
    %Get the ratio of baseline to the data for further inference
    ratio_baseline_to_data{dt} = (baseline_to_be_added{dt}*(size(edges,2)-1))/size(reconstit_arr{dt},1);
    
    %Now, get the 'partial' offset, which is used for fitting. This is
    %caleld 'TruthoffsetPartial', as this is based on data known 'before'
    %the fitting via MLE.
    %Probably, both of these have something to do with the mean bin size,
    %rather than these strange formulae correcting for maxdist...
    if linorlogBGsubtract == 'lin'
        truthoffsetpartial{dt} = (2e-6/maxdist)*(baseline_to_be_added{dt}*bgbinningnr)/...
            (size(reconstit_arr{dt},1)-(baseline_to_be_added{dt}*bgbinningnr));
    elseif linorlogBGsubtract == 'log'
        bins = logspace(log10(minlogpoint),log10(maxdist),bgbinningnr);
        truthoffsetpartial{dt} = (2e-6/maxdist)*(baseline_to_be_added{dt}*bgbinningnr)/...
            (size(reconstit_arr{dt},1)-(baseline_to_be_added{dt}*bgbinningnr))*mean(bins);
        clear bins
    end
end
%Calculate the number of entries for the JD array --> This should be
%equal to roughly what's in the true input - based on the number of JDs
%found and the BG ratio.
for dt = 1:size(bgratios,2)
    nrentries_JDarray(dt) = ceil(size(JDarrSignalCell{dt},1)*(1-bgratios(dt)));
    JDarraydata{dt} = reconstit_arr{dt}(randperm(size(reconstit_arr{dt},1),nrentries_JDarray(dt)));
end
JDonlydata = {JDarraydata;reconstit_arr};