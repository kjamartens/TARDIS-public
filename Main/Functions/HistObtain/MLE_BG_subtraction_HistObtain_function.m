function [histInfo,time,BGcurve_interp,bestratio,JDarrBG,JDarrSignalCell,signalCurve_interp,stepsizearraytrue,BGarrtotsize] = MLE_BG_subtraction_HistObtain_function(poslist,frame_dist_BG,maxdist,frame_dist_Hist_arr,bgbinningnr,startpointBG,linorlogBGsubtract,minlogpoint,visualisation,verbose,callfromUI)
%% Function to obtain information for histograms based on TARDIS
% Required inputs:
% poslist           1-x-4 Array with frame-x-y-z position (positions in meters)
% frame_dist_BG     Frame distance (array) used to find background value. Normally used values: 0 or larger than maximum length of track
% max_dist          Maximum distance used for fitting (in meters)
% frame_time        Frame time duration in seconds
% frame_dist_Hist_arr         Arry of dt bins (i.e. [1 2 3])
% verbose           True for verbose output, False for silent operation
%---------------------------------------------------------------------------------------------------------------------------------------
% Obtained outputs
% histInfo          Obtained histogram information after BG subtraction
% time              Duration of function in seconds
% BGcurve_interp    Information of the background curve
%%
if verbose
    dispUIorCommandWindow(' ',callfromUI);
    dispUIorCommandWindow('-----------------------------------------',callfromUI);
    dispUIorCommandWindow(['Starting HistObtain_TARDIS ' linorlogBGsubtract],callfromUI);
end
tic
%For debug:
%correct_linkage_calculation(poslist,frame_dist_Hist_arr,poslistunordered)
%% Variable initiation
mindist_onlynoise = startpointBG;
maxdist_onlynoise = maxdist;

%% Obtaining BG information
if verbose; dispUIorCommandWindow('Starting BG obtain',callfromUI); end
%Doing multiple frame_dist_BG if wanted
BGarrtotsize = [];
for k = 1:size(frame_dist_BG,2)
    %Terrible slow implementation atm, but who cares
    BGarr_maxdist_singledt{k} = JD_Rel(poslist,frame_dist_BG(k),maxdist,0); %Provides the Background relative JD array
    BGarr_maxdist_singledt{k}(BGarr_maxdist_singledt{k}==0) = []; %Remove zeros from the BG array (only happens if frame_dist_BG = 0, and the localization JD is calculated with respect to themselves)
    % NOTE: I THINK THIS LINE CAN BE SAFELY REMOVED! BGarr_maxdist = BGarr(BGarr<=maxdist); %limit BGarr to the maximum distance specified
    BGarrtotsize(k) = size(BGarr_maxdist_singledt{k},1);
end
%Get all BGarr_maxdist_singledt in a single BGarr_maxdist
BGarr_maxdist = zeros(sum(BGarrtotsize),1);
for k = 1:size(frame_dist_BG,2)
    startpoint = 1;
    if k > 1
        startpoint = sum(BGarrtotsize(1:k-1))+1;
    end
    endpoint = startpoint+BGarrtotsize(k)-1;
    BGarr_maxdist(startpoint:endpoint) = BGarr_maxdist_singledt{k};
end

%% Obtain 'signal' JDs
if verbose; dispUIorCommandWindow('Starting Signal obtain',callfromUI); end
%Obtain the jumping distances (JD) for all combinations based on the dt
%values
inputarr = []; %iniate an empty array
tempsteparr = cell(1,size(frame_dist_Hist_arr,1));
size_dt = zeros(size(frame_dist_Hist_arr,1),1);
for dt = frame_dist_Hist_arr %Loop over all dt bins wanted
    tempsteparr{dt} = JD_Rel(poslist(:,1:3),dt,maxdist,0); %Get the relative JDs for every dt bin
    inputarr = [inputarr; tempsteparr{dt}]; %Put in a large array combining all dt bins for later use
    size_dt(dt) = size(tempsteparr{dt},1); %Get the number of JD values for every dt bin
    
    %Next few lines are cleanup, but should normally never happen
    if size(tempsteparr{dt}(tempsteparr{dt}<0),1) > 0
        tempsteparr{dt}(tempsteparr{dt}<0) = []; %Cleanup, remove negative values
        disp('Something strange happening with cleanup');
    end
    if size(tempsteparr{dt}(tempsteparr{dt}>maxdist),1) > 0
        tempsteparr{dt}(tempsteparr{dt}>maxdist) = []; %Cleanup, remove too large v alues
        disp('Something strange happening with cleanup');
    end
    %Cleanup of zeros
    tempsteparr{dt}(tempsteparr{dt}==0) = [];
end
%% Get information about BG to put in
if verbose; dispUIorCommandWindow('Starting Interpolation',callfromUI); end
% The background information is created as follows: A certain number of
% interpolation bins are created, linearly spanning 0 to maxdist.
% Essentially, a histogram of the BG values is then created on these
% interpolation bins, and the values are stored. Then, for all signal
% values, the corresponding BG value is found by getting a closest match of
% this BG histogram (x-position), and interpolation between the closest
% positions (linearly).

% The result (output_BG_alldt) is a value for how likely it is for all
% positions to belong to the BG rather than to the signal. This could be
% visaulised by for example scatter(inputarr,output_BG_alldt). Note that it
% not rescaled or normalised.
%Note that there are 2 options here, based on log or lin fitting
if linorlogBGsubtract == 'lin'
    if ~isempty(BGarr_maxdist)
        %BGcurve_interp gives a x,y list of interpolated points
        a = linspace(0,maxdist,bgbinningnr+1);
        get_points = zeros(1,size(a,2)-1);
        for i = 1:size(a,2)-1
            get_points(i) = a(i)*0.5+a(i+1)*0.5;
        end
        BGcurve_interp = interpolate_BGCurve_arr(bgbinningnr,get_points,BGarr_maxdist);
    else
        dispUIorCommandWindow('Empty background array! All tracks are completely separated temporally',callfromUI);
        %Make BG arrays linear, i.e. unaffecting the data
        BGcurve_interp(:,1) = linspace(0,maxdist,figurebinnr+1);
        BGcurve_interp(:,2) = ones(size(BGcurve_interp,1),1);
    end
elseif linorlogBGsubtract == 'log'
    if ~isempty(BGarr_maxdist)
        %BGcurve_interp gives a x,y list of interpolated points
        a = logspace(log10(minlogpoint),log10(maxdist),bgbinningnr+1);
        get_points = zeros(1,size(a,2)-1);
        for i = 1:size(a,2)-1
            get_points(i) = sqrt(a(i)*a(i+1));
        end
        BGcurve_interp = interpolate_BGCurve_arr_log(bgbinningnr,get_points,BGarr_maxdist,minlogpoint,maxdist);
    else
        dispUIorCommandWindow('Empty background array! All tracks are completely separated temporally',callfromUI);
        %Make BG arrays linear, i.e. unaffecting the data
        BGcurve_interp(:,1) = logspace(log10(minlogpoint),log10(maxdist),figurebinnr+1);
        BGcurve_interp(:,2) = ones(size(BGcurve_interp,1),1);
    end
else
    dispUIorCommandWindow(['Error with variable linorlogBGsubtract! Should be "lin" or "log", but it is: ' linorlogBGsubtract],callfromUI);
end
%output_BG_alldt is a list for all inputarr values on how likely it is that
%it belongs to the BG. This is later rescaled according to a fraction of
%the BG inside the MLE
% output_BG_alldt = pdfBGFunction(inputarr,BGcurve_interp);

%% Do same as BG interpolation, but for signal
signalcurve_interp = cell(1,size(frame_dist_Hist_arr,1));
for dt = frame_dist_Hist_arr
    if ~isempty(tempsteparr{dt})
        %BGcurve_interp gives a x,y list of interpolated points
        %Note that get_points is on purpose 'recycled' from BG get_points
        if linorlogBGsubtract == 'lin'
            signalcurve_interp{dt} = interpolate_BGCurve_arr(bgbinningnr,get_points,tempsteparr{dt});
        elseif linorlogBGsubtract == 'log'
            signalcurve_interp{dt} = interpolate_BGCurve_arr_log(bgbinningnr,get_points,tempsteparr{dt},minlogpoint,maxdist);
        end
    end
end

%% Find ratio based on last x positions
clear bestratio
%pre-allocate array
bestratio = zeros(size(frame_dist_Hist_arr,1),1);
%Loop over dt bins
for dt = frame_dist_Hist_arr
    clear t rmse ratiocounter
    %Have a very over-sampled ratiocounter to test all of them. Script is
    %fast enough
    ratiocounter = [0:0.001:1];
    %pre-allocate (and reset) RMSE-array
    rmse = zeros(size(ratiocounter,2),1);
    %also calculating what ratio has the mean at zero, but unused (don't
    %think it helps, very small deviations anyhow
    rmseZeroMean = zeros(size(ratiocounter,2),1);
    %Only select the xdata positions between the mindist and maxdist
    %specified
    xdatarelevantmin = find(BGcurve_interp(:,1)>=mindist_onlynoise);
    xdatarelevantmax = find(BGcurve_interp(:,1)<=maxdist_onlynoise);
    %Combine min and max
    xdatarelevant = [];
    for i = 1:size(xdatarelevantmin,1)
        if ~isempty(find(xdatarelevantmin(i) == xdatarelevantmax,1))
            xdatarelevant = [xdatarelevant; xdatarelevantmin(i)];
        end
    end
    %Test all ratiovals, and calculate the RMSE
%%
% figure(5);clf(5);
% subplot(1,2,1)
% plot(BGcurve_interp(:,1),BGcurve_interp(:,2),'k-x')
% set(gca,'XScale','log')
% subplot(1,2,2)
% plot(signalcurve_interp{dt}(:,1),signalcurve_interp{dt}(:,2),'k-x')
% set(gca,'XScale','log')
%%
    for ratioval = 1:size(ratiocounter,2)
        t = signalcurve_interp{dt}(xdatarelevant,2) - BGcurve_interp(xdatarelevant,2)*ratiocounter(ratioval);
        rmse(ratioval) = sqrt(sum(t.^2));
        
        t = signalcurve_interp{dt}(xdatarelevant,2) - BGcurve_interp(xdatarelevant,2)*ratiocounter(ratioval);
        rmseZeroMean(ratioval) = sum(t);
    end
    %Find the minimum RMSE for this dt bin and save it
    [~,b] = min(rmse);
    bestratio(dt) = ratiocounter(b);
    [~,b] = min(abs(rmseZeroMean));
    bestratioZeroMean(dt) = ratiocounter(b);
end
if verbose; dispUIorCommandWindow(['BG ratio found at values ' num2str(bestratio)],callfromUI); end
%ZeroMean unused, but saved in case. Time is neglicable I assume
% if verbose; disp(['BG ratioZeroMean found at values ' num2str(bestratioZeroMean)]); end
%% Get true JD dist
%Error clause to only do this if the information is given (normally not, so
%NORMALLY an empty stepsizearraytrue matrix is given)
if exist('poslistunordered')
    stepsizearraytrue = cell(1,size(frame_dist_Hist_arr,1));
    for dt = frame_dist_Hist_arr
        list = [1:9:45000];
        poslistundered_forstepsize = poslistunordered(:,[1 2 3]);
        for i = 1:5000
            poslistundered_forstepsize(list(i):list(i)+8,4) = i;
        end
        stepsizearraytrue{dt} = tracks2steps(poslistundered_forstepsize,dt);
        stepsizearraytrue{dt}(stepsizearraytrue{dt}==0) = [];
    end
else
    for dt = frame_dist_Hist_arr
        stepsizearraytrue{dt}=[];
    end
end

%% Create output
%Create empty cell
histInfo = cell(1,size(frame_dist_Hist_arr,1));
clear t
%Loop over dts
for dt = frame_dist_Hist_arr
    %Subtract BG from signal based on ratio found
    t = signalcurve_interp{dt}(:,2) - BGcurve_interp(:,2)*bestratio(dt);
    %Normalize the output and store it in histInfo - which is the output of
    %the script
    t = t./sum(t);
    histInfo{dt} = [BGcurve_interp(:,1) t];
end
    %%
    %Write some variables for output
    JDarrBG = BGarr_maxdist;
    JDarrSignalCell = tempsteparr;
    signalCurve_interp = signalcurve_interp;
time = toc;
if verbose; dispUIorCommandWindow(['Ending HistObain_TARDIS, total time: ' num2str(time) 's'],callfromUI); dispUIorCommandWindow('-----------------------------------------',callfromUI); end
end
%% Functions and stuff
function [output] = interpolate_BGCurve_arr(bgbinningnr,bins,BGlist)
%Function to create an interpolated BG curve. The idea is to get a
%bgbinningnr-by-2 list of x,y positions of the BG curve. Therefore, the
%BGlist is simply histogrammed, and the histogram positions are further
%interpolated makima-style.
%Get histogram from BGlist
[h,edges] = histcounts(BGlist,bgbinningnr,'Normalization','probability');
%Get the middle of the bins
midbins = zeros(1,size(edges,2)-1);
for i = 1:size(midbins,2)
    midbins(i) = (edges(i)+edges(i+1))/2;
end
%Interpolate histogram from BGlist to be used later
try
    output(:,2) = interp1([0,midbins],[0,h],bins,'makima');%spline is ok
    output(:,1) = bins;
catch
    keyboard
end
end


function [output] = interpolate_BGCurve_arr_log(bgbinningnr,bins,BGlist,bgmindist,bgmaxdist)
%Function to create an interpolated BG curve. The idea is to get a
%bgbinningnr-by-2 list of x,y positions of the BG curve. Therefore, the
%BGlist is simply histogrammed, and the histogram positions are further
%interpolated makima-style.
%Get histogram from BGlist
[h,edges] = histcounts(BGlist,logspace(log10(bgmindist),log10(bgmaxdist),bgbinningnr),'Normalization','probability');
%Get the middle of the bins
midbins = zeros(1,size(edges,2)-1);
for i = 1:size(midbins,2)
    midbins(i) = sqrt(edges(i)*edges(i+1));
end
%Interpolate histogram from BGlist to be used later
try
    output(:,2) = interp1([0,midbins],[0,h],bins,'makima');%spline is ok
    output(:,1) = bins;
catch
    keyboard
end
end
%%
function stepsizearray = tracks2steps(tracklist,framediff)
%Expected input: track with [frame, x, y, trackID]
stepsizearray = zeros(length(tracklist)-length(unique(tracklist(:,4))),1);
counter = 1;
%create accumarray
acarr = accumarray(tracklist(:,4),1);
%loop over all tracks
trackids = unique(tracklist(:,4));
for n = 1:length(trackids)
    i = trackids(n);
    startid = 1;
    if i > 1
        startid = sum(acarr(1:i-1))+1;
    end
    endid = length(tracklist);
    if i < length(acarr)
        endid = sum(acarr(1:i));
    end
    for k = startid+framediff:endid
        stepsizearray(counter) = sqrt((tracklist(k-framediff,2)-tracklist(k,2))^2+(tracklist(k-framediff,3)-tracklist(k,3))^2);
        counter = counter+1;
    end
end
end
%% Calculate wrong/correct linkages - only used for debugging - requires ,poslistunordered
function theorratio = correct_linkage_calculation(poslist,frame_dist_Hist_arr,poslistunordered)
ae = accumarray(poslist(:,1),1);
totallinks = zeros(size(ae,1),size(frame_dist_Hist_arr,1));
for dt = frame_dist_Hist_arr
    for fr1 = 1:size(ae,1)-dt
        yposstart = sum(ae(1:fr1-1))+1;
        posfr1 = poslist(yposstart:yposstart+ae(fr1)-1,:);
        posfr2 = poslist(yposstart+ae(fr1):yposstart+ae(fr1)+ae(fr1+dt)-1,:);
        totallinks(fr1,dt) = size(posfr1,1)*size(posfr2,1);
    end
end
%Now for correct linkages

clear stepsizearraytrue
stepsizearraytrue = cell(1,size(frame_dist_Hist_arr,1));
corrlinks = zeros(size(frame_dist_Hist_arr,1),1);
for dt = frame_dist_Hist_arr
    list = [1:9:45000];
    poslistundered_forstepsize = poslistunordered(:,[1 2 3]);
    for i = 1:5000
        poslistundered_forstepsize(list(i):list(i)+8,4) = i;
    end
    stepsizearraytrue{dt} = tracks2steps(poslistundered_forstepsize,dt);
    stepsizearraytrue{dt}(stepsizearraytrue{dt}==0) = [];
    corrlinks(dt) = size(stepsizearraytrue{dt},1);
end

theorratio = 1-(corrlinks)./sum(totallinks,1);
end