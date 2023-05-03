%Expected input:
%poslist: frame-x-y-z list (m)
%frame_dist_BG:
function [paramEsts, tottime, time, bgarr] = URDA_HO_MLE_function(poslist,settingsURDA)
tic
%% Obtain variables from settings and warn user if some are not provided
%Initial variables
frame_dist_BG = 9; % Which dt's are used to get BG (f.e. [9:1:20])
startpointBG = 1e-6; %what JD the BG-only area is starting
maxdist = 2e-6; %What the max JD is
bgbinningnr = 100; %nr of BG bins
dt_arr = [1 2 3]; %which dt's are used to get signal
vis = 1; %Visibility of BG subtracting results
verbose = 1; 
fitvis = 1; %Visibility of fit results
linorlogBGsubtract = 'lin'; %'lin' or 'log' - later probably also 'dda'
linorlogvis = linorlogBGsubtract;
minlogpoint = 10^(-8.5); %starting point if using log BG subtract
frame_time = 0.01; %Frame time in seconds
loc_unc = 15e-9; %Localization uncertainty in meters
norm_bins = 10000; %Nr of bins used for normalization - everything > 1000 seems fine
start = [5e-12]; %Startpoint of fitting
lb = [1e-16]; %LB of fit
ub = [30e-12]; %UB of fit
debug = false; %Debug info

%String list of all variables
variablearr = {'frame_dist_BG','startpointBG','maxdist','bgbinningnr','dt_arr',...
    'vis','verbose','fitvis','linorlogBGsubtract','linorlogvis','minlogpoint','frame_time','loc_unc','norm_bins',...
    'start','lb','ub','debug'};

%Check for every variable if it's in the settingsURDA struct. If not, an
%warning is given and the pre-set value is used.
for i = 1:size(variablearr,2)
    if ~isfield(settingsURDA,variablearr{i})
        disp(['No settingsURDA.' variablearr{i} ' provided! Using a pre-set value of ' num2str(eval(variablearr{i}))]);
    else
        eval([variablearr{i} ' = settingsURDA.' variablearr{i} ';']);
    end
end

%% Obtain Histogram from JD
[output,time,bgarr,bgratios] = MLE_BG_subtraction_HistObtain_function(poslist,frame_dist_BG,maxdist,dt_arr,bgbinningnr,startpointBG,linorlogBGsubtract,minlogpoint,vis,verbose);%,poslistunordered);

%% Create 'fake' JD array
%Get number of entries in total
totalreconstit_arr_entries = 10000; % This is a sort of minimum? It's slightly higher than this because of removal of negatives
%Enable or disable having a 'baseline' for every bin. If this is disabled,
%results are likely biased
baseline_per_bin = true;
endpoint = maxdist;
%Find Edges of the histogram
if linorlogBGsubtract == 'lin'
    startpoint = 0;
    halfwidthsize = (output{1}(2,1)-output{1}(1,1))/2;
    if (halfwidthsize-output{1}(1,1)) > 1e-16
        disp('This is going to fail! Probably change code')
    end
    %For linear, simply interpolate between start and endpoint
    edges = [startpoint:halfwidthsize*2:endpoint];
elseif linorlogBGsubtract == 'log'
    %For log, interpolate in log-space
    loghalfwidthsize = (log10(output{1}(2,1))-log10(output{1}(1,1)))/2;
    edges = [log10(minlogpoint):loghalfwidthsize*2:log10(endpoint)];
    edges = 10.^edges;
end
%Loop over the dt bins
for dt = dt_arr
    %Create previous and next interpoints, which dictate the shape of how
    %the random points in the reconstit_arr should be arranged.
    %F.e. a bin that's part of a quickly rising bin, shouldn't be a
    %'block', but more like a slope.
    %Idea is to get nr of JDs in the reconstit_arr from the output, but the
    %shape from prev and next-interp.
    for i = 1:size(edges,2)-1
        if i == 1
            previnterp(i) = 0;
        elseif i == size(edges,2)-1
            nextinterp(i) = 0;
        end
        
        if i > 1
            previnterp(i) = 0.5*(output{dt}(i,2)+output{dt}(i-1,2));
        end
        if i < size(edges,2)-1
            nextinterp(i) = 0.5*(output{dt}(i,2)+output{dt}(i+1,2));
        end
    end
    %Figure for debugging:
    if debug
        if linorlogBGsubtract=='lin'
            figure(90);clf(90);
            hold on
            plot(output{dt}(:,1),output{dt}(:,2),'k-x');
            plot(output{dt}(:,1)-halfwidthsize,previnterp,'r-o');
            plot(output{dt}(:,1)+halfwidthsize,nextinterp,'r-x');
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
            nrentries_bin_check(i) = floor(output{dt}(i,2)*totalreconstit_arr_entries);
        end
        %find the baseline that should be added - should be a positive nr 1
        %higher than the smallest negative bin in the histogram
        baseline_to_be_added{dt} = max(0,min(nrentries_bin_check)*-1+1);
    end
    for i = 1:size(edges,2)-1
        %For every bin, get the number of entries based on the size of the
        %bin and the baseline
        if ~baseline_per_bin
            nrentries_bin = max(round(output{dt}(i,2)*totalreconstit_arr_entries),0);
        else
            nrentries_bin = round(output{dt}(i,2)*totalreconstit_arr_entries)+baseline_to_be_added{dt};
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
            fancysubbins = 1;
            if fancysubbins
                %Create the nr of subbins
                subbins = 100;
                %Get a probability of a localization being in the subbins
                subbinprob = linspace(previnterp(i),nextinterp(i),subbins);
                subbinprob = subbinprob./sum(subbinprob);
                subbinprobcum = cumsum(subbinprob);
                tarr = zeros(nrentries_bin,1);
                randarr = rand(nrentries_bin,1);
                %Loop over all entries, choose a subbin based on the
                %cumulative probabilities, and throw a random value in this
                %subbin
                for cntr = 1:nrentries_bin
                    chosenbin = find((subbinprobcum-randarr(cntr))>0,1,'first');
                    tarr(cntr) = edges(i)+((halfwidthsize*2)/subbins)*(chosenbin-1)+rand(1,1)*(halfwidthsize*2)/subbins;
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
    % Visualisation of blockyness due to binsize for debugging
    if debug
        if linorlogBGsubtract=='lin'
            figure(91);clf(91);
            mf = 4;
            ratio_baseline_to_data{dt} = (baseline_to_be_added{dt}*(size(edges,2)-1))/size(reconstit_arr{dt},1);
            histogram(reconstit_arr{1},linspace(0,maxdist,mf*(bgbinningnr+1)),'Normalization','probability')
            hold on
            yd = output{1}(:,2)*(1-ratio_baseline_to_data{1})+ones(size(output{1}(:,2)))./bgbinningnr*ratio_baseline_to_data{1};
            plot(output{1}(:,1),(yd./sum(yd))/mf,'r-x');
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

%% Start the actual MLE fit
dtbinsize = [];
reconstit_arr_MLE=[];
%Add the data of all dt's behind each other
for dt = dt_arr
    dtbinsize = [dtbinsize size(reconstit_arr{dt},1)];
    reconstit_arr_MLE = [reconstit_arr_MLE; reconstit_arr{dt}];
end
%Get a numeric log/lin variable
switch linorlogBGsubtract
    case 'lin'
        numericvaluelinorlog = 1;
    case 'log'
        numericvaluelinorlog = 2;
end
%Initiate the MLE fit
custompdf = @(x,D)MLEfitter(x,D,truthoffsetpartial,loc_unc,frame_time,dtbinsize,norm_bins,numericvaluelinorlog,debug,minlogpoint);
mleoptions = statset('mlecustom');
mleoptions.MaxIter = 2000; %Maximum iterations
mleoptions.MaxFunEvals = mleoptions.MaxIter; %Maximum function evaluations
if verbose; disp('Starting MLE fit'); end
%Performing the MLE fit
paramEsts = mle(reconstit_arr_MLE, 'pdf',custompdf, 'start',start, ...
    'lower',lb, 'upper',ub, 'options',mleoptions);
if verbose; disp(['MLE fit completed with params: ' num2str(paramEsts*1e12)]); end

%% Visualisation fit
if fitvis
    figure(3);clf(3);
    %Loop over dt bins
    for dt = dt_arr
        %Some info to check baseline - unused
%         xdata = [0:1e-8:maxdist];
%         frame_timedt = frame_time*dt;
%         MLEt = 1e-12+loc_unc^2/frame_timedt;
%         yd1 = xdata.*exp(-(xdata.^2)./(4*MLEt*frame_timedt));
%         yd1 = yd1./sum(yd1);
%         yd2 = ones(size(yd1));
%         yd2 = yd2./sum(yd2);
%         ydata = yd1.*(1-ratio_baseline_to_data{dt})+yd2.*(ratio_baseline_to_data{dt});
%         ydata = ydata./sum(ydata);
%         %     figure(6);clf(6);
%         %     plot(xdata,ydata);
%         %     hold on
%         foundoffsetblack{dt} = ydata(end).*(size(xdata,2)/100);
%         DMLE = paramEsts(1);
%         DMLE = 1e-12;
%         DMLE=DMLE+loc_unc^2/frame_time;
%         ydata = xdata.*exp(-(xdata.^2)./(4*DMLE*frame_time))+DMLE*1e4*truthoffsetpartial{dt}*dt;
%         foundoffsetcyan{dt} = (DMLE*1e4*truthoffsetpartial{dt}*dt)./sum(ydata).*(size(xdata,2)/100);
        
        % Actual visualization
        subplot(1,max(dt_arr),(dt));
        %Create histogram and xdata for fit plot
        if linorlogBGsubtract == 'lin'
            h=histogram(reconstit_arr{dt},linspace(0,maxdist,bgbinningnr+1),'Normalization','Probability');
            xdata = [0:1e-8:maxdist];
        elseif linorlogBGsubtract == 'log'
            h=histogram(reconstit_arr{dt},logspace(log10(minlogpoint),log10(maxdist),bgbinningnr+1),'Normalization','Probability');
            xdata = logspace(log10(minlogpoint),log10(maxdist),bgbinningnr+1);
        end
        histFormat(h);
        hold on
        %Get some variables for the fitted curve, such as frame_time for
        %this dt bin, and loc_unc-corrected DiffCoeff
        frame_timedt = frame_time*dt;
        DMLE = paramEsts(1);
        DMLE=DMLE+loc_unc^2/frame_timedt;
        %Different visalization for lin or log
        if linorlogBGsubtract == 'lin'
            offset = DMLE*1e4*truthoffsetpartial{dt}*dt;
            ydata = xdata.*exp(-(xdata.^2)./(4*DMLE*frame_timedt))+offset;
%             foundoffsetcyan{dt} = (DMLE*1e4*truthoffsetpartial{dt}*dt)./sum(ydata).*(size(xdata,2)/bgbinningnr);
            plot(xdata,ydata./sum(ydata).*(size(xdata,2)/bgbinningnr),'k-','LineWidth',2,'DisplayName',['MLE: ' num2str(paramEsts*1e12)]);
        elseif linorlogBGsubtract == 'log'
            offset = (DMLE*1e4*dt)*truthoffsetpartial{dt};
            ydata = (xdata.*exp(-(xdata.^2)./(4*DMLE*frame_timedt))).*xdata+offset;
            plot(xdata,(ydata./sum(ydata).*(size(xdata,2)/bgbinningnr)),'k-','LineWidth',2)
        end
        %Get maxheight to later scale the subplots
        maxheight(dt) = max(max(ydata),max(h.Values));
        %Plot the baseline
        yd2 = output{dt}(:,2)*totalreconstit_arr_entries+baseline_to_be_added{dt};
        plot(output{dt}(:,1),ones(size(output{dt}(:,1)))./sum(ones(size(output{dt}(:,1))))*(baseline_to_be_added{dt}/mean(yd2)),'k--')
        %Add a title
        title(['Data minus BG with fit, dt = ' num2str(dt)])
    end
    %Rescale the axis and do something nice with axes notation
    for dt = dt_arr
        subplot(1,max(dt_arr),dt)
        if dt == min(dt_arr)
            ylabel('Occurance')
        else
            set(gca,'YTick',[]);
        end
        xlabel('JD (m)')
        
        if linorlogBGsubtract == 'lin'
            axis([0 maxdist 0 max(maxheight)*1.1])
        elseif linorlogBGsubtract == 'log'
            set(gca,'XScale','log')
            axis([minlogpoint maxdist 0 max(maxheight)*1.1])
        end
    end
end
tottime = toc;
end
%%
function [output] = MLEfitter(xdata,D,offsetpartial,loc_unc,frame_time,dtbinsize,norm_bins,linorlogBGsubtract,debug,minlogpoint)
%assuming that f.e. dtbins = [1 2 3], meaning 1,2,3 frame delay, and all
%the same x-range/xdata
%Output will be [ydatadt1 ydatadt2 ydatadt3]

%Pre-allocate output
output = zeros(size(xdata));
%Loop over dt
for dt = 1:size(dtbinsize,2)
    %Adept frame time to reflect dt
    frame_time_dt = frame_time*dt;
    %calculate offset on the fly. Unsure yet about exact formula
    if linorlogBGsubtract == 1
        offset{dt} = (D+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
    elseif linorlogBGsubtract == 2 %ABSOLUTELY NO CLUE WHY
        offset{dt} = (D+loc_unc^2/frame_time)*1e4*offsetpartial{dt}*dt;
    end
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
        outputsingledt{dt} = pdfSinglePopFunction(xdatadt,D,offset{dt},loc_unc,frame_time_dt);
        %Normalize for area
        outputsingledt{dt} = normalize_area(xdatadt,outputsingledt{dt},D,offset{dt},loc_unc,frame_time_dt,norm_bins);
    elseif linorlogBGsubtract == 2 %Logarithmic
%         keyboard
        outputsingledt{dt} = pdfSinglePopFunction_log(xdatadt,D,offset{dt},loc_unc,frame_time_dt);
        %Normalize for area
        outputsingledt{dt} = normalize_area_log(xdatadt,outputsingledt{dt},D,offset{dt},loc_unc,frame_time_dt,norm_bins);
    end
    %Avoid some errors (never called in good fittings)
    outputsingledt{dt}(outputsingledt{dt}<=0) = 1e-150;
    
    output(startpointxdata:endpointxdata) = outputsingledt{dt};
end

%DEBUG drawing - can later be expanded/cleaned in some nice stuff
if debug
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
            Dplot = D+loc_unc^2/frame_time_dt;
            offsetplot = offset{dt};
            xdatapl = linspace(0,max(xdata),100);
            ydata = (xdatapl.*exp(-(xdatapl.^2)./(4*Dplot*frame_time_dt)))+offsetplot;
            plot(xdatapl,(ydata./sum(ydata).*(size(xdatapl,2)/100)),'k-','LineWidth',2,'DisplayName',num2str(D))
        elseif linorlogBGsubtract == 2 %Logarithmic
            h=histogram(xdata(startpointxdata:endpointxdata),logspace(log10(minlogpoint),log10(max(xdata)),100),'Normalization','probability');
            histFormat(h);
            set(gca,'XScale','log')
            Dplot = D+loc_unc^2/frame_time_dt;
            offsetplot = offset{dt};
            xdatapl = logspace(log10(minlogpoint),log10(max(xdata)),100);
            ydatal = (xdatapl.*exp(-(xdatapl.^2)./(4*Dplot*frame_time_dt))).*xdatapl+offsetplot;
            plot(xdatapl,(ydatal./sum(ydatal).*(size(xdatapl,2)/100)),'k-','LineWidth',2,'DisplayName',num2str(D))
        end
        xlabel('JD (m)')
        ylabel('Occurance')
        legend()
        title(['DEBUG fitting in progress - dt:' num2str(dt)])
    end
    drawnow
end
end
%% Normalize area scripts - lin and log
%No real difference detected in lin or log normalization...
function [output] = normalize_area(xdata,input,D,offset,loc_unc,frame_time,norm_bins)
%Calculate the area by pseudo-integration. Calculate the area for
%many small bins and correct output for the sum of those bins.
%In principle, the area of the curve specified by the parameters (D,
%loc_unc, f) is calculated, and the input is normalized for the sum of this
%area
area_calc_bin_size = (max(xdata)-min(xdata))/norm_bins; %Make 1000 bins
check_range = min(xdata):area_calc_bin_size:(max(xdata)-area_calc_bin_size); %Specifiy x-range it should be checked for
mean_x_pos_bin = (check_range+area_calc_bin_size/2)'; %Get middle x position of every bin
%Get area of all bins via the same function as in the regular probability
%calculation. This gives the total area of the function
%Different function for BG or population
mean_y_pos_bin = pdfSinglePopFunction(mean_x_pos_bin,D,offset,loc_unc,frame_time);
%Get total area
area_bin = mean_y_pos_bin*area_calc_bin_size;
totarea = sum(area_bin);
%Correct for total area
output = input./totarea;
end

function [output] = normalize_area_log(xdata,input,D,offset,loc_unc,frame_time,norm_bins)
%Calculate the area by pseudo-integration. Calculate the area for
%many small bins and correct output for the sum of those bins.
%In principle, the area of the curve specified by the parameters (D,
%loc_unc, f) is calculated, and the input is normalized for the sum of this
%area
separataion_calc_bin_size_log = (log10(max(xdata))-log10(min(xdata)))/norm_bins; %Make 1000 bins
check_range_log = log10(min(xdata)):separataion_calc_bin_size_log:(log10(max(xdata))-separataion_calc_bin_size_log); %Specifiy x-range it should be checked for
mean_x_pos_bin_logtolin = 10.^((check_range_log+separataion_calc_bin_size_log/2)'); %Get middle x position of every bin

logspacing = logspace(log10(min(xdata)),log10(max(xdata)),norm_bins+1);
area_calc_bins_log = logspacing(2:end)-logspacing(1:end-1);
%Get area of all bins via the same function as in the regular probability
%calculation. This gives the total area of the function
%Different function for BG or population
mean_y_pos_bin = pdfSinglePopFunction_log(mean_x_pos_bin_logtolin,D,offset,loc_unc,frame_time);
%Get total area
area_bin = mean_y_pos_bin.*area_calc_bins_log';
totarea = sum(area_bin);
%Correct for total area
output = input./totarea;
end
%% Actual fitting description - lin and log variants
function output = pdfSinglePopFunction(xdata,D,offset,loc_unc,frame_time)
output = xdata.*exp(-(xdata.^2)./(4*(D+loc_unc^2/frame_time)*frame_time))+offset;
end

function output = pdfSinglePopFunction_log(xdata,D,offset,loc_unc,frame_time)
D=D+loc_unc^2/frame_time;
output =  xdata.*exp(-(xdata.^2)./(4*D*frame_time))+offset./xdata;
end
%% Formatting function
function histFormat(h)
h.EdgeColor = 'k';
h.EdgeAlpha = 0.2;
h.FaceColor = [0.5 0.5 0.5];
h.FaceAlpha = 0.6;
end
%% Backup for scripts to go to die
% function MLE - unknown offset variant, unused for now
% function [output] = MLEfitter_unknownoffset(xdata,D,offset,loc_unc,frame_time,dtbinsize,norm_bins,linorlogBGsubtract)
% %Pre-allocate output
% output = zeros(size(xdata));
% %Loop over dt
% for dt = 1:size(dtbinsize,2)
%     %Adept frame time to reflect dt
%     frame_time_dt = frame_time*dt;
%     if dt > 1
%         startpointxdata = sum(dtbinsize(1:dt-1))+1;
%     else
%         startpointxdata = 1;
%     end
%     endpointxdata = startpointxdata+dtbinsize(dt)-1;
%     xdatadt = xdata(startpointxdata:endpointxdata);
%     %Perform fitting on this dt-separated data
%     if linorlogBGsubtract == 1 %Linear
%         outputsingledt{dt} = pdfSinglePopFunction(xdatadt,D,offset,loc_unc,frame_time_dt);
%         %Normalize for area
%         outputsingledt{dt} = normalize_area(xdatadt,outputsingledt{dt},D,offset,loc_unc,frame_time_dt,norm_bins);
%     elseif linorlogBGsubtract == 2 %Logarithmic
%         outputsingledt{dt} = pdfSinglePopFunction(xdatadt,D,offset,loc_unc,frame_time_dt);
%         %Normalize for area
%         outputsingledt{dt} = normalize_area_log(xdatadt,outputsingledt{dt},D,offset,loc_unc,frame_time_dt,norm_bins);
%     end
%     %Avoid some errors (never called in good fittings)
%     outputsingledt{dt}(outputsingledt{dt}<=0) = 1e-150;
%     
%     %Create output-array
%     output(startpointxdata:endpointxdata) = outputsingledt{dt};
% end
% end

% LSQ fit
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% frame_time = 0.01;
% loc_unc = 15e-9*1e6;
% dtbins = dt_arr;
% fun = @(D,xdata)LSQfitter(xdata,D,loc_unc,frame_time,dt_arr);
%
% startD = [2e-12]*1e12;
% clear options
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','OptimalityTolerance',1e-12,'FunctionTolerance',1e-9);%,'StepTolerance',1e-15);
% % options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','OptimalityTolerance',1e-12);%,'FunctionTolerance',1e-9);%,'StepTolerance',1e-15);
% lb = [0.01e-12]*1e12;
% ub = [10e-12]*1e12;
%
% fitxdata = repmat(output{1}(:,1)*1e6,size(dt_arr,2),1);
% fitydata = zeros(size(fitxdata));
% for dt = 1:size(dt_arr,2)
%     startp = 1+size(output{1}(:,2),1)*(dt-1);
%     endp = startp+size(output{1}(:,2),1)-1;
%     fitydata(startp:endp) = output{dt}(:,2);
% end
% clear startp endp
% outp = lsqcurvefit(fun,startD,fitxdata,fitydata,lb,ub,options);
% toc
% outp
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% function LSQ - UNUSED
% function [output] = LSQfitter(xdata,D,loc_unc,frame_time,dt_arr)
% %assuming that f.e. dt_arr = [1 2 3], meaning 1,2,3 frame delay, and all
% %same x-range/xdata
% %Output should be [ydatadt1 ydatadt2 ydatadt3]
% % keyboard
% xdatasingle = xdata(1:size(xdata,1)/size(dt_arr,2));
%
% for dt = dt_arr
%     psuedoframetime = dt.*frame_time;
%     outputsingledt{dt} = xdatasingle.*exp(-(xdatasingle.^2)./(4*(D+loc_unc^2/psuedoframetime)*psuedoframetime));
%     outputsingledt{dt} = outputsingledt{dt}./sum(outputsingledt{dt});
% end
%
% output = zeros(size(xdata));
% for dt = 1:size(dt_arr,2)
%     startp = 1+size(outputsingledt{1},1)*(dt-1)
%     endp = startp+size(outputsingledt{1},1)-1
%     output(startp:endp) = outputsingledt{dt};
% end
%
% end


% Poor-mans fitting
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% clear output_Diff
% Darr = [0.1e-12:0.001e-12:2e-12];
% for Dcntr = 1:size(Darr,2)
%     D = Darr(Dcntr);
%     loc_unc = 15e-9;
%     frame_time = 0.01;
%     D=D+loc_unc^2/frame_time;
%     xdata = output{1}(:,1);
%     output_Diff(Dcntr,:) =  xdata.*exp(-(xdata.^2)./(4*D*frame_time));%exp(-(xdata./(4*D*frame_time)));
%     output_Diff(Dcntr,:) = output_Diff(Dcntr,:)./sum(output_Diff(Dcntr,:));
%     rmse(Dcntr) = sum(sqrt((output_Diff(Dcntr,:)'-output{1}(:,2)).^2));
% end
% [a,b] = min(rmse);
% minDfound = Darr(b)
% output_Diff = output_Diff./(size(output{1},2)/size(output_Diff,2));
