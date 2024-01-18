%Function that combines HistObtain (HO) method to get initial fit
%parameters, and uses those for a full fit afterwards.
%Koen J.A. Martens et al., 2021.

%Expected input:
%poslist: frame-x-y-z list (m)
function [parameters, parametersCI, paramEsts, HOfitCI, tottime, time, bgarr, truthoffsetpartial, Visualisation_HO_outputCell,Visualisation_FF_outputCell,anaDDAvisInfoHO,anaDDAvisInfoFF,JDonlydata,SwiftParameters,MHTparameters] = URDA_HO_FF_function(poslist,settingsTARDIS)
%% Check for required toolboxes
%Get the toolbox info
MATLABinfo = ver;
%Find whether the required toolboxes are installed, and throw an error
%otherwise
if ~any(strcmp({MATLABinfo.Name}, 'Curve Fitting Toolbox'))
    warning('Curve Fitting Toolbox not found!');
    error('Curve Fitting Toolbox not installed! Install it via Home - Add-Ons, or check mathworks.com/products/curvefitting.html')
end
if ~any(strcmp({MATLABinfo.Name}, 'Signal Processing Toolbox'))
    warning('Signal Processing Toolbox not found!');
    error('Signal Processing Toolbox not installed! Install it via Home - Add-Ons, or check mathworks.com/products/signal.html')
end
if ~any(strcmp({MATLABinfo.Name}, 'Optimization Toolbox'))
    warning('Optimization Toolbox not found!');
    error('Optimization Toolbox not installed! Install it via Home - Add-Ons, or check mathworks.com/products/optimization.html')
end
if ~any(strcmp({MATLABinfo.Name}, 'Statistics and Machine Learning Toolbox'))
    warning('Statistics and Machine Learning Toolbox not found!');
    error('Statistics and Machine Learning Toolbox not installed! Install it via Home - Add-Ons, or check mathworks.com/products/statistics.html')
end
clear MATLABinfo
%%
tic;
%Sort poslist on frame
poslist = sortrows(poslist,1);
%% Obtain variables from settings and warn user if some are not provided
%Initial variables
frame_time = 0.01; %- %Frame time in seconds
loc_unc = 15e-9; %- %Localization uncertainty in meters
frame_dist_BG = [20:40]; %- % Which dt's are used to get BG (f.e. [9:1:20])
startpointBG = 1e-6; %- %what JD the BG-only area is starting
maxdist = 3e-6; %- %What the max JD is
bgbinningnr = 100; %- %nr of BG bins
linorlogBGsubtract = 'log'; %'lin' or 'log'
minlogpoint = 10^(-8.5); %starting point if using log BG subtract - used to be -8.5
dt_arr = [1 2 3]; %- %which dt's are used to get signal
debug = 0; %Debug info
vis = debug; %Visibility of BG subtracting results
verbose = 1;
verboseReal = verbose; %actual verbose of fitting etc
fitvisHO = debug; %Visibility of fit results
linorlogvis = linorlogBGsubtract; %- ;
norm_bins = 50000; %Nr of bins used for normalization - everything > 1000 seems fine
start_1pop = '[0.1-20]'; %- %Startpoint of fitting between 1e-11 and 1e-13
lb_1pop = '[1e-5]'; %LB of fit
ub_1pop = '[50]'; %UB of fit
start_2pop = '[0.1-20 0.1-20 0.1-10]';%-  %Startpoint of fitting between 1e-11 and 1e-13
lb_2pop = '[1e-5 1e-5 1e-5]'; %LB of fit
ub_2pop = '[50 50 100]'; %UB of fit
start_3pop = '[0.05 1 29 0.9]'; %Startpoint of fitting between 1e-11 and 1e-13
lb_3pop = '[1e-5 1e-5 1e-5 0]'; %LB of fit
ub_3pop = '[1e4 1e-4 1e4 100]'; %UB of fit
start_aDDA = '[10-80 10-80 0.1-10]'; %- %Startpoint of anaDDAfit
lb_aDDA = '[1e-5 1e-5 1e-5]';
ub_aDDA = '[1e4 1e4 1e2]';
visualisationMLEIntFit = 0; %Visualisation of MLE fit at the end
performestimationfit = 1; %- %Perform estimation fit
performsecondfit = 1; %Perform second (better) MLE fit
populations = 1; %- %number of populations to be fit. Make this 0 to fit anaDDA!
callfromUI = 0; %Set to the UI app (callfromUI = app) if called from UI - usefull for text output;
createJDsonly = 0; %Only create a JD list rather than try to fit something
fitWithBleach = 0; %fit BG/populations with bleach kinetics rather than with arbitrary values
fixRatios2pop = 1; %Fit 2pop with fixed ratios (true) or not (false)
AutoChoosePop = 0; %Automatically choose between 1, 2, or 3 (not yet implemented 3) populations
stroboFrameTime = 0; %Stroboscopic frame time. Set to -1 for full-frame, or to a value in seconds.
freefit_locunc = 0; %EXPERIMENTAL! Only implemented in log_BG subtracting. Set to 1 to fit the loc unc rather than provide. If this is 1, then the provided loc unc will be used as initial guess
StoreSWIFTparameters = 0; %Store SWIFT JSON as output. Needs to be called from UI for correct folder storage.
StoreMHTparameters = 0; %Store MHT (multiple-hypothesis-tracking) .xml as output. Needs to be called from UI for correct folder storage.
noiseDensity = 0; %Value that is required for extracting swift params. Should be the density of the noise in loc/um2/frame
customPDF = 0; %Information about a custom PDF, or 0 if none.
customPDF_start = []; %Start, upper, lower bounds of a custom PDF
customPDF_upper = [];
customPDF_lower = [];
verboseMLEIntFit = verboseReal;
AlternativeLookupPosList = 0; %Have a N-by-4 alternative lookup list (similar width as poslist) if you want to have the 'lookup' values different than the 'start' values. Make 0 if this isn't applicable.

%String list of all variables
variablearr = {'frame_dist_BG','startpointBG','maxdist','bgbinningnr','dt_arr',...
    'vis','verbose','fitvisHO','linorlogBGsubtract','linorlogvis','minlogpoint','frame_time','loc_unc','norm_bins',...
    'debug','verboseMLEIntFit','verboseReal','visualisationMLEIntFit','performestimationfit','performsecondfit','populations',...
    'start_1pop','lb_1pop','ub_1pop','start_2pop','lb_2pop','ub_2pop','start_3pop','lb_3pop','ub_3pop',...
    'start_aDDA','callfromUI','createJDsonly','fitWithBleach','fixRatios2pop','AutoChoosePop','stroboFrameTime','freefit_locunc',...
    'StoreSWIFTparameters','StoreMHTparameters','noiseDensity','customPDF','customPDF_start','customPDF_upper','customPDF_lower','AlternativeLookupPosList'};

%Check for every variable if it's in the settingsTARDIS struct. If not, an
%warning is given and the pre-set value is used.
for i = 1:size(variablearr,2)
    if ~isfield(settingsTARDIS,variablearr{i})
        disp(['No settingsTARDIS.' variablearr{i} ' provided! Using a pre-set value of ' num2str(eval(variablearr{i}))]);
        try
            eval(['settingsTARDIS.' variablearr{i} '=' num2str(eval(variablearr{i})) ';']);
        catch
            eval(['settingsTARDIS.' variablearr{i} '=[' num2str(eval(variablearr{i})) '];']);
        end
    else
        eval([variablearr{i} ' = settingsTARDIS.' variablearr{i} ';']);
    end
end

% %% EXPERIMENTAL SHOULD BE REMOVED
% dt_arr = [3,4];
%     dispUIorCommandWindow('Starting estimation fit',callfromUI);


%Empty output because it might crash
anaDDAvisInfoHO = [];
anaDDAvisInfoFF = [];
SwiftParameters = '';
MHTparameters = '';

%Make callfromUI empty if zero inputted for later easy-code
if callfromUI == 0
    callfromUI = [];
end

[~, neworder] = sort(lower(fieldnames(settingsTARDIS)));
settingsTARDIS = orderfields(settingsTARDIS, neworder);

%% Setting pdfSettings
pdfSettings.frame_time = frame_time;
pdfSettings.loc_unc = loc_unc;
pdfSettings.populations = populations;
%Note: size_dt following later
pdfSettings.verbose = verboseMLEIntFit;
pdfSettings.fitWithBleach = fitWithBleach;
pdfSettings.fixRatios2pop = fixRatios2pop;
pdfSettings.freefit_locunc = freefit_locunc;
if stroboFrameTime == -1
    stroboFrameTime = frame_time;
end
pdfSettings.strobo_frame_time = stroboFrameTime;
%% Deal with random start positions
% Note that even though MATLAB thinks the start values are unused, they are
% called for in an eval() function.
[settingsTARDIS, start_1pop, start_2pop, start_aDDA] = RandomStartPosToNormal_UI(settingsTARDIS);

%Make upper/lower limits values rather than text
lb_1pop = eval(settingsTARDIS.lb_1pop);
ub_1pop = eval(settingsTARDIS.ub_1pop);
lb_2pop = eval(settingsTARDIS.lb_2pop);
ub_2pop = eval(settingsTARDIS.ub_2pop);
lb_aDDA = eval(settingsTARDIS.lb_adda);
ub_aDDA = eval(settingsTARDIS.ub_adda);
%% Obtain Histogram from JD
% Gain information on the JD arrays
% frame_dist_BG=0
try
    [histInfo,time,bgarr,bgratios,JDarrBG,JDarrSignalCell,signalCurve_interp,stepsizearraytrue,BGarrtotsize] = ...
        MLE_BG_subtraction_HistObtain_function(poslist,frame_dist_BG,maxdist,dt_arr,bgbinningnr,startpointBG,linorlogBGsubtract,minlogpoint,0,verbose,callfromUI,'AlternativeLookupPosList',AlternativeLookupPosList);%,poslistunordered);
catch e
    %Print error message, Issue1 bugfix
    fprintf(1,'There was an error! The message was:\n%s\n',e.message);
    %     keyboard
end
% keyboard
if vis
    additionalinfo.Data.settingsURDA = settingsTARDIS;
    additionalinfo.Data.anaDDAvisInfoHO.size_dt_HO = settingsTARDIS.dt_arr;
    %     keyboard
    Visualisation_HO_BGsubtract('','',bgbinningnr,dt_arr,JDarrBG,maxdist,minlogpoint,linorlogBGsubtract,bgarr,JDarrSignalCell,...
        signalCurve_interp,stepsizearraytrue,bgratios,startpointBG,additionalinfo);
    clear additionalinfo
end
Visualisation_HO_outputCell.bgbinningnr = bgbinningnr;
Visualisation_HO_outputCell.dt_arr = dt_arr;
Visualisation_HO_outputCell.JDarrBG = JDarrBG;
Visualisation_HO_outputCell.maxdist = maxdist;
Visualisation_HO_outputCell.minlogpoint = minlogpoint;
Visualisation_HO_outputCell.linorlogBGsubtract = linorlogBGsubtract;
Visualisation_HO_outputCell.bgarr = bgarr;
Visualisation_HO_outputCell.JDarrSignalCell = JDarrSignalCell;
Visualisation_HO_outputCell.signalCurve_interp = signalCurve_interp;
Visualisation_HO_outputCell.stepsizearraytrue = stepsizearraytrue;
Visualisation_HO_outputCell.bgratios = bgratios;
Visualisation_HO_outputCell.startpointBG = startpointBG;
Visualisation_HO_outputCell.histInfo = histInfo;
%Testing this
for k = 1:size(histInfo,2)
    histInfo{k}(histInfo{k}(:,2)<0,2) = 0;
    histInfo{k}(:,2) = histInfo{k}(:,2)./(sum(histInfo{k}(:,2)));
    %     dispUIorCommandWindow('histInfo contained < 0 values - correcting for this!',callfromUI)
end
%% Create 'reconstituted' JD array
[JDonlydata, reconstit_arr, truthoffsetpartial, baseline_to_be_added] = ...
    createReconstitutedJDarray(JDarrSignalCell,bgratios,linorlogBGsubtract,histInfo,debug,maxdist,minlogpoint,dt_arr,bgbinningnr);
%%

% plot([1:3],2.^(-[1:3]/3))
%Get SWIFT parameters via corrected BG ratios
%Also saves the SWIFT JSON if requested via UI
% % % SwiftParameters = getCorrectedBGRatio(JDonlydata,JDarrBG,reconstit_arr,JDarrSignalCell,frame_dist_BG,poslist,settingsTARDIS);
%% Swift and MHT options
%Check if we want SWIFT output, because we can approximate those
%try to get the bleach half-time here from bgratios - not perfect
fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[4,max(1-bgratios)]);
ft = fittype('b*2^(-(x+1)/a)','options',fo);
fitc = fit(dt_arr',[1-bgratios]',ft);
bleach_time_half = (fitc.a-.5)*settingsTARDIS.frame_time;

dispUIorCommandWindow(['Approximate bleach half-time: ' num2str(bleach_time_half)],callfromUI);

parameters = bleach_time_half; %This is the only parameter that obtain_swift_params looks for :)
if settingsTARDIS.StoreMHTparameters || settingsTARDIS.StoreSWIFTparameters
    [SwiftParameters,MHTparameters] = obtain_swift_params(parameters,JDarrSignalCell,poslist,maxdist,frame_dist_BG,settingsTARDIS,noiseDensity,BGarrtotsize,JDonlydata,callfromUI);

    if settingsTARDIS.StoreSWIFTparameters
        dispUIorCommandWindow(['Found blink ratio for swift: ' num2str(SwiftParameters.p_blink)],callfromUI);
    end
end
%% JD-only creation if wanted
%If we only want the JD list, we have all the information we want, and we
%can stop the analysis here.
if createJDsonly
    %Return variables so function doesn't crash
    anaDDAvisInfoHO=[];anaDDAvisInfoFF=[];
    time = toc; parameters = [];parametersCI=[];paramEsts=[];HOfitCI=[];tottime = time;Visualisation_FF_outputCell=[];
    return
end

%% Start the actual HO_MLE fit
% This fits the BG-subtracted JD histograms, and assumes that there is no
% BG present. This is a nice way to get initial values.
% Also known as 'estimation fit' or similar

%Firstly, not done if autochoosepop is used.
if AutoChoosePop
    performestimationfit = 0;
    dispUIorCommandWindow('Estimation fit not done because of AutoChoosePop!',callfromUI)
end

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
%Initiate the HO fit
if verbose && performestimationfit
    dispUIorCommandWindow('Starting estimation fit',callfromUI);
end
%Set MLE options
mleoptions = TARDIS_setMLEoptions('HO');
%%
% if ~useAIC %If AIC is not to be used
if populations > 0 %If populations are meant to be fitted rather than anaDDA
    %Now perform normal fit with info from AIC, or from beginning
    %         pdf_HO_pops = @(x,varargin)HistOnly_MLEFit_estimation(x,truthoffsetpartial,populations,loc_unc,frame_time,dtbinsize,norm_bins,numericvaluelinorlog,debug,minlogpoint,verbose,callfromUI,varargin);
    pdf_HO_pops = @(x,varargin)HistOnly_MLEFit_estimation(x,truthoffsetpartial,pdfSettings,dtbinsize,norm_bins,numericvaluelinorlog,debug,minlogpoint,callfromUI,varargin);
    %Performing the MLE fit
    start = eval(eval(['start_' num2str(populations) 'pop']));
    eval(['lb = lb_' num2str(populations) 'pop;']);
    eval(['ub = ub_' num2str(populations) 'pop;']);
    if performestimationfit %don't do it if there's no esimation fit
        %remove zeros from reconstit_arr and set them to a low number
        reconstit_arr_MLE(reconstit_arr_MLE==0) = 1e-15;
        [paramEsts, HOfitCI] = mle(reconstit_arr_MLE, 'pdf',pdf_HO_pops, 'start',start, ...
            'lower',lb, 'upper',ub, 'options',mleoptions);
    else
        paramEsts = []; HOfitCI = [];
    end
elseif populations == 0 %Else fit anaDDA to the data
    %Import information for anaDDA2.0
    %         filename_temp = fullfile('C:\Users\Koen Martens\Documents\Martens_WorkingDirectory\Scientific\MATLAB\TARDIS_testingADDA20', 'locdisttable.mat');
    try
        %When running from MATLAB, do this...
        load('.\Main\Functions/anaDDA_requiredScripts_V20/locdisttable.mat'); %Loads pre-generated correlated localization error distributions
    catch
        %Otherwise, when installed as exe, do this...
        tempA = load('locdisttable.mat','locdist');
        locdist = tempA.locdist;
        tempB = load('locdisttable.mat','rangex'); %Loads pre-generated correlated localization error distributions
        rangex = tempB.rangex;
    end
    settingsTARDIS.inputAnaDDA.rangex = rangex; %rangex loaded line above
    for k = 2:8
        settingsTARDIS.inputAnaDDA.locdist{k} = griddedInterpolant(rangex(:,k),locdist(:,k));
    end
    %Get some input over from orig settings
    settingsTARDIS.inputAnaDDA.sigmaerror = settingsTARDIS.loc_unc*1e6; %in um
    settingsTARDIS.inputAnaDDA.trackingwindow = settingsTARDIS.maxdist;%*1e3; %careful: hard-coded number
    settingsTARDIS.inputAnaDDA.frametimerange = settingsTARDIS.frame_time.*[1:size(settingsTARDIS.dt_arr,2)];

    %End import
    anaDDAvisInfoHO = [];
    anaDDAvisInfoFF = [];

    %Get size_dt
    [size_dt,~] = getSizesInputFromJDCell(JDarrSignalCell);
    pdfSettings.size_dt = size_dt;

    %Quite some info now is also used for the FF fit later...
    try
        inputAnaDDA = settingsTARDIS.inputAnaDDA;
%         keyboard
        %Testing random operation on some input to test if all is
        %loaded
        %Testing a few things 2023.01.31:
%         inputAnaDDA.compensatetracking = 1;
        inputAnaDDA.Steptime = 1e-5;
        inputAnaDDA.upperDfree*2;
        disp('Loaded anaDDAsettings correctly')
    catch
        inputAnaDDA = Generateinputfile_KMmodified;
        inputAnaDDA.numberofspecies = 1;
        inputAnaDDA.upperDfree = 10; % in um2/s
        inputAnaDDA.sigmaerror = 15e-3; %in um
        inputAnaDDA.confinement = 0; %for now
        inputAnaDDA.radiusofcell = 3000; %careful: hard-coded number
        inputAnaDDA.lengthcell = 3000; %careful: hard-coded number
        inputAnaDDA.compensatetracking = 0;
        inputAnaDDA.trackingwindow = 3000; %careful: hard-coded number
        inputAnaDDA.fixedparameters(:) = -1;
        inputAnaDDA.framerange = 1; %hardcoded and should always be this way
        inputAnaDDA.precision = 2^16; %Normally 50000
        disp('Loaded anaDDAsettings from script')
    end

    %Always get this one - shouldn't be necessarily user-defined, and
    %is important to get right:
    inputAnaDDA.frametimerange = frame_time.*[1:size(dt_arr,2)];
    inputAnaDDA.fitlocerror = 0;
    %Next line(s) testing out anaDDA v2.0
    %         inputAnaDDA.fixedparameters = [1 -1 -1 -1 0; -1 -1 -1 -1 0;-1 -1 -1 -1 0];  % Fixed parameters for fitting. Each row is a species, the columns are fraction, koff, kon, Dfree and D1. -1 if not fixed(Default: [1 -1 -1 -1 0; -1 -1 -1 -1 0;-1 -1 -1 -1 0];)
    inputAnaDDA.integrationinterval = 2000; %changed from 200 31.01
    %         load('C:\Users\Koen Martens\Documents\GitHub\URDA\anaDDA_requiredScripts_V2\locdisttable.mat');
    %         inputAnaDDA.locdist = locdist;
    %         inputAnaDDA.rangex = rangex;
    %D has to be an array of Diff Values - ones - frame time. Frame time should
    %scale with dt values
    D_alldtbins = [];
    rangeD_alldtbins = [];

    inputarr_D = [];
    for dt = 1:size(dt_arr,2)
        %Doing something weird with loc unc
        %             keyboard
        tempsteparr_D{dt} = (JDarrSignalCell{dt}.*1e6).^2./(4*frame_time*dt);
%                     tempsteparr_D{dt} = ((JDarrSignalCell{dt}.*1e6).^2 + (loc_unc*1e6)^2/(dt))./(4*frame_time*dt);
        inputarr_D = [inputarr_D; tempsteparr_D{dt}]; %Put in a large array combining all dt bins for later use
    end

    for dt = 1:size(dt_arr,2)
        D_singledtbin{dt} = [tempsteparr_D{dt} ones(size(tempsteparr_D{dt},1),1).*dt ones(size(tempsteparr_D{dt},1),1).*frame_time.*dt];

        %value for D max - alllmost constant over all dt bins
        maxDfree_singledtbin{dt} = inputAnaDDA.upperDfree+inputAnaDDA.sigmaerror^2/min(frame_time.*dt);
        maxrangeD_singledtbin{dt} =-log(maxDfree_singledtbin{dt}*1e-10)*maxDfree_singledtbin{dt};
        %rangeD then appears to be a log range from a low value to maxDfree?
        rangeD_singledtbin{dt} =maxrangeD_singledtbin{dt}/(inputAnaDDA.precision*2):maxrangeD_singledtbin{dt}/inputAnaDDA.precision:maxrangeD_singledtbin{dt};

        %Aggregate matrices in one big matrix
        D_alldtbins = [D_alldtbins;D_singledtbin{dt}];
        rangeD_alldtbins = [rangeD_alldtbins;rangeD_singledtbin{dt}];

        %Pre-calculate locerrorpdf and locerrorpdfcorrected. THese should be
        %passed on to GeneratefixedDdistributions_KM and DDistributiongenerator

        %Now some PDFs are calculated, which scale with frametimes.
        locerror_singledtbin{dt} = inputAnaDDA.sigmaerror.^2;%/(frame_time.*dt);%
        inputAnaDDA.frametime = (frame_time.*dt);
        [locerrorpdf_singledtbin{dt},locerrorpdfcorrected_singledtbin{dt}] = makelocerrordistributions(rangeD_singledtbin{dt},locerror_singledtbin{dt},inputAnaDDA);

        % Create determined parameters Dfixed, fitspecies, and fixedspecies
        [Dfixed_singledtbin{dt}, fitspecies_singledtbin{dt}, fixedspecies_singledtbin{dt}] = GeneratefixedDdistributions_KM(inputAnaDDA.numberofspecies, inputAnaDDA.fixedparameters, rangeD_singledtbin{dt},locerrorpdfcorrected_singledtbin{dt},inputAnaDDA);
    end
    inputAnaDDA.frametime = (frame_time.*1); %reset - added 31.01
    %Using nomenclature from anaDDA for easyness
    Dlistdata = D_alldtbins';
    numberofspecies = inputAnaDDA.numberofspecies;
    fixedparameters = inputAnaDDA.fixedparameters;

    % meanD of dataset
    meanD = mean(Dlistdata(1,:));
    fixedparameters(1,1) = 1;
    % Number of degrees of freedom (all the ones that are 0)
    % depend on species and which parameters are fixed
    fixedparameterstemp = fixedparameters(1:numberofspecies,:);

    %Setting parameters
    % startparameters = [0.3 50 30 8;0.3 1 1 1; 0.3 1 1 1];
    indexfittingparameters = fixedparameterstemp==-1;
    indexfittingparameters = [indexfittingparameters; zeros(3-numberofspecies,4)];
    indexfittingparameters = logical(indexfittingparameters);

    % startparam = [startparameters(indexfittingparameters); 0.5];
    % lowerbound = [zeros(length(startparam)-1,1)+0.000001;0];
    % upperbound = [ones(sum(fixedparameterstemp(:,1)==-1),1) ;input.upperstartkon*input.upperstartkoff*ones(sum(sum(indexfittingparameters(:,2:3)==1)),1);10*ones(sum(indexfittingparameters(:,4)==1),1); 1];

    %Rephrasing parameters
    Frametimelist = Dlistdata(3,:);
    Numberofframes = Dlistdata(2,:);
    Dlistdata = Dlistdata(1,:);
%     doing some strange sorting here, unsure why.
    [Numberofframes, sortind]=sort(Numberofframes);
    Dlistdata = Dlistdata(sortind);
    Frametimelist = Frametimelist(sortind);
    [Frametimelist, sortind]=sort(Frametimelist);
    Dlistdata = Dlistdata(sortind);
    Numberofframes = Numberofframes(sortind);
    %         keyboard
    for j = 1:numel(inputAnaDDA.frametimerange)
        %     table = tabulate(Numberofframes(Frametimelist==input.frametimerange(j)));
        frequency(j,:) = size_dt(j);%table(:,2);
        inputAnaDDA.frametime = inputAnaDDA.frametimerange(j);
        [fx(:,j),fy(:,j)] = Generateconfinedfunction(0:0.05:inputAnaDDA.upperDfree,rangeD_alldtbins,inputAnaDDA);
    end
    %Reset to 1 frametime delay
    inputAnaDDA.frametime = inputAnaDDA.frametimerange(1);

    for dt = 1:size(dt_arr,2)
        if inputAnaDDA.compensatetracking == true
            maxD = (inputAnaDDA.trackingwindow*inputAnaDDA.pixelsize)^2/(4*inputAnaDDA.frametime);
            maxDindtracking{dt} = (maxD./(rangeD_alldtbins(dt,2)-rangeD_alldtbins(dt,1)));
        else
            maxDindtracking{dt} = 0;
        end
    end

    % Get BG curve in D
    if ~isempty(JDarrBG)
        for dt = 1:size(size_dt,2)
            JDarrBG(JDarrBG==0) = []; %Remove zeros from the BG array (only happens if frame_dist_BG = 0, and the localization JD is calculated with respect to themselves)
            %                 BGarr_maxdist_D{dt} = ((JDarrBG).*1e6).^2./(4*frame_time*dt);
            BGarr_maxdist_D{dt} = ((JDarrBG.*1e6).^2)./(4*frame_time*dt);
            maxdist_D{dt} = (maxdist.*1e6).^2./(4*frame_time*dt);
            logminval_D{dt} = (minlogpoint.*1e6).^2./(4*frame_time*dt);
            bgbinningnr_interpolate = settingsTARDIS.bgbinningnr;
            BGcurve_interp_D{dt} = interpolate_BGCurve_log(maxdist_D{dt},logminval_D{dt},bgbinningnr,bgbinningnr_interpolate,BGarr_maxdist_D{dt});
        end
    else
        disp('Empty background array! All tracks are completely separated temporally')
        %Make BG arrays linear, i.e. unaffecting
        maxdist_D = (maxdist.*1e6).^2./(4*frame_time*dt);
        logminval_D = (minlogpoint.*1e6).^2./(4*frame_time*dt);
        BGcurve_interp_D(:,1) = logspace(log10(logminval_D),log10(maxdist_D),bgbinningnr+1);
        BGcurve_interp_D(:,2) = ones(size(BGcurve_interp_D,1),1);
        %Set lower upper start bounds to account for no bg
        %     startparam,lowerbound,
        %Probably to fix?
        %     startparam(end-size(size_dt,2)+1:end) = 0.0005;
        %     lowerbound(end-size(size_dt,2)+1:end) = 0;
        %     upperbound(end-size(size_dt,2)+1:end) = 0.001;
    end
    %         keyboard
    for dt = 1:size(size_dt,2)
        output_BG_alldt_D{dt} = pdfBGFunction_log(inputarr_D,BGcurve_interp_D{dt});
    end

    %Everything above here is also important for the FF fit later!!

    %Now perform an anaDDA fit with linear background on the
    %reconstructed histogram.
    %First get D info from the reconstit JD arr
    inputReconstitarr_D = [];
    for dt = 1:size(dt_arr,2)
        tempsteparrReconstit_D{dt} = (reconstit_arr{dt}.*1e6).^2./(4*frame_time*dt);
        inputReconstitarr_D = [inputReconstitarr_D; tempsteparrReconstit_D{dt}]; %Put in a large array combining all dt bins for later use
    end
    %Then get a linear BG datalist
    for dt = 1:size(size_dt,2)
        output_BG_alldt_D_linear{dt} = pdfBGFunction_log(inputReconstitarr_D,BGcurve_interp_D{dt});
        %             output_BG_alldt_D_linear{dt} = ones(size(output_BG_alldt_D_linear{dt}));
        BGcurve_interp_D_linear{dt}(:,1) = BGcurve_interp_D{dt}(:,1);
        BGcurve_interp_D_linear{dt}(:,2) = ones(size(BGcurve_interp_D{dt}(:,2)));
        BGcurve_interp_D_linear{dt}(:,2) = BGcurve_interp_D_linear{dt}(:,2)./sum(BGcurve_interp_D_linear{dt}(:,2));
    end

    D_alldtbins_reconstit=[];
    for dt = 1:size(dt_arr,2)
        D_singledtbin_reconstit{dt} = [tempsteparrReconstit_D{dt} ones(size(tempsteparrReconstit_D{dt},1),1).*dt ones(size(tempsteparrReconstit_D{dt},1),1).*frame_time.*dt];

        %Aggregate matrices in one big matrix
        D_alldtbins_reconstit = [D_alldtbins_reconstit;D_singledtbin_reconstit{dt}];
    end

    %Using nomenclature from anaDDA for easyness
    Dlistdata_reconstit = D_alldtbins_reconstit';

    % meanD of dataset
    meanD_reconstit = mean(Dlistdata_reconstit(1,:));

    %Rephrasing parameters
    Frametimelist = Dlistdata_reconstit(3,:);
    Numberofframes_reconstit = Dlistdata_reconstit(2,:);
    Dlistdata_reconstit = Dlistdata_reconstit(1,:);
    %doing some strange sorting here, unsure why.
    [Numberofframes_reconstit, sortind]=sort(Numberofframes_reconstit);
    Dlistdata_reconstit = Dlistdata_reconstit(sortind);
    Frametimelist = Frametimelist(sortind);
    [Frametimelist, sortind]=sort(Frametimelist);
    Dlistdata_reconstit = Dlistdata_reconstit(sortind);
    Numberofframes_reconstit = Numberofframes_reconstit(sortind);


    for dt = 1:size(reconstit_arr,2) %Loop over all dt bins wanted
        size_dt_reconstit(dt) = size(reconstit_arr{dt},1);
    end
    frequency_reconstit = size_dt_reconstit';

    %Now perform a similar fit as done in the FF fit later
    custompdfReconstitAnaDDA = @(Dlistdata_reconstit,varargin) pdfAnaDDAMLE_multidt(inputReconstitarr_D,BGcurve_interp_D_linear, Numberofframes_reconstit,Frametimelist,...
        inputAnaDDA,rangeD_alldtbins,fitspecies_singledtbin{1}, fixedspecies_singledtbin{1}, ...
        (((Dfixed_singledtbin{dt}(1)))), fixedparameters,indexfittingparameters,fx,fy,maxDindtracking,...
        frequency_reconstit,size_dt_reconstit,(locerrorpdfcorrected_singledtbin),output_BG_alldt_D_linear,bgbinningnr,fitWithBleach,verbose,callfromUI,varargin);

    mleoptions = TARDIS_setMLEoptions('HO');

    %         upperbound = [500 500 2.01 ones(1,max(dt_arr))*1]';
    %

    if performestimationfit
        startparam = [eval(eval('start_aDDA')) ones(1,max(dt_arr))*0.005]';
        lowerbound = [lb_aDDA ones(1,max(dt_arr))*0]';
        upperbound = [ub_aDDA ones(1,max(dt_arr))*0.01]';
        %             startparam(1:3) = rand(3,1).*[70 70 3]'+[10 10 1]'; %10-80, 10-80, 1-4
        %                 keyboard
        %aDDA estimation fit
        [paramEsts] = mle(Dlistdata_reconstit, 'pdf',custompdfReconstitAnaDDA,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound,'Options',mleoptions);
        paramEsts
        HOfitCI=[];

        %Visualisation info for HO-anaDDA
        [~,anaDDAvisInfoHO.visual] = pdfAnaDDAMLE_multidt(inputReconstitarr_D,BGcurve_interp_D_linear, Numberofframes_reconstit,Frametimelist,...
            inputAnaDDA,rangeD_alldtbins,fitspecies_singledtbin{1}, fixedspecies_singledtbin{1}, ...
            (((Dfixed_singledtbin{dt}(1)))), fixedparameters,indexfittingparameters,fx,fy,maxDindtracking,...
            frequency_reconstit,size_dt_reconstit,(locerrorpdfcorrected_singledtbin),output_BG_alldt_D_linear,bgbinningnr,fitWithBleach,verbose,callfromUI,paramEsts);
        anaDDAvisInfoHO.DlistHO = inputReconstitarr_D;
        anaDDAvisInfoHO.size_dt_HO = size_dt_reconstit;
    end
end
% end
if performestimationfit
    if populations == 1
        if verbose
            stringval = sprintf('First fit (estimate) completed with DiffCoeff %.2f (%.2f-%.2f).\n',paramEsts(1),HOfitCI(1,1),HOfitCI(2,1));
            dispUIorCommandWindow(stringval,callfromUI);
            dispUIorCommandWindow(' ',callfromUI);
        end
    end
    if populations == 2
        if verbose
            stringval = sprintf('First fit (estimate) completed with DiffCoeffs: %.2f (%.2f-%.2f) and %.2f (%.2f-%.2f) at ratio 1:%.2f (%.2f-%.2f).\n',paramEsts(1),HOfitCI(1,1),HOfitCI(2,1),paramEsts(2),HOfitCI(1,2),HOfitCI(2,2),paramEsts(3),HOfitCI(1,3),HOfitCI(2,3));
            dispUIorCommandWindow(stringval,callfromUI);
            dispUIorCommandWindow(' ',callfromUI);
        end
    end

    %% Visualisation HO part
    if fitvisHO
        figure(4);clf(4);
        %Loop over dt bins
        for dt = dt_arr
            totalreconstit_arr_entries = size(reconstit_arr{dt},1);
            %Some info to check baseline - unused
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
            %Different visalization for lin or log
            if populations > 0 %Not doing this for aDDA
                if linorlogBGsubtract == 'lin'
                    if populations == 1
                        DMLE = paramEsts(1)/1e12+loc_unc^2/frame_timedt;
                        offset = DMLE*1e4*truthoffsetpartial{dt}*dt;
                        ydata = xdata.*exp(-(xdata.^2)./(4*DMLE*frame_timedt))+offset;
                        %             foundoffsetcyan{dt} = (DMLE*1e4*truthoffsetpartial{dt}*dt)./sum(ydata).*(size(xdata,2)/bgbinningnr);
                        plot(xdata,ydata./sum(ydata).*(size(xdata,2)/bgbinningnr)*Visualisation_HO_outputCell.bgratios(1)/Visualisation_HO_outputCell.bgratios(dt),'k-','LineWidth',2,'DisplayName',['MLE: ' num2str(paramEsts*1e12)]);
                        %Get maxheight to later scale the subplots
                        maxheight(dt) = max(max(ydata),max(h.Values));
                    elseif populations == 2
                        %slightly different: plot both on top of baseline
                        %get baseline
                        yd2 = histInfo{dt}(:,2)*totalreconstit_arr_entries+baseline_to_be_added{dt};
                        baseline = ((baseline_to_be_added{dt}/mean(yd2))/size(histInfo{dt}(:,1),1))*ones(size(xdata));

                        DMLE1 = paramEsts(1)/1e12+loc_unc^2/frame_timedt;
                        offset1 = DMLE1*1e4*truthoffsetpartial{dt}*dt;
                        ydata1 = xdata.*exp(-(xdata.^2)./(4*DMLE1*frame_timedt))+offset1;
                        ydata1norm = (ydata1./sum(ydata1).*(size(xdata,2)/bgbinningnr))*(1/(1+paramEsts(3)));
                        ydata1norm = ydata1norm-ydata1norm(end)+baseline;
                        plot(xdata,ydata1norm*Visualisation_HO_outputCell.bgratios(1)/Visualisation_HO_outputCell.bgratios(dt),'b-','LineWidth',1,'DisplayName',['MLE: ' num2str(paramEsts*1e12)]);

                        DMLE2 = paramEsts(2)/1e12+loc_unc^2/frame_timedt;
                        offset2 = DMLE2*1e4*truthoffsetpartial{dt}*dt;
                        ydata2 = xdata.*exp(-(xdata.^2)./(4*DMLE2*frame_timedt))+offset2;
                        ydata2norm = (ydata2./sum(ydata2).*(size(xdata,2)/bgbinningnr))*(paramEsts(3)/(1+paramEsts(3)));
                        ydata2norm = ydata2norm-ydata2norm(end)+baseline;
                        plot(xdata,ydata2norm*Visualisation_HO_outputCell.bgratios(1)/Visualisation_HO_outputCell.bgratios(dt),'b-','LineWidth',1,'DisplayName',['MLE: ' num2str(paramEsts(1)*1e12)]);

                        ydatasum = ydata1norm-baseline+ydata2norm-baseline + baseline;
                        plot(xdata,ydatasum*Visualisation_HO_outputCell.bgratios(1)/Visualisation_HO_outputCell.bgratios(dt),'k-','LineWidth',1,'DisplayName',['MLE: ' num2str(paramEsts(2)*1e12) '-' num2str(paramEsts(3))]);

                        %Get maxheight to later scale the subplots
                        maxheight(dt) = max(max(ydatasum),max(h.Values));
                    end
                elseif linorlogBGsubtract == 'log'
                    if populations == 1
                        DMLE = paramEsts(1)/1e12+loc_unc^2/frame_timedt;
                        offset = (DMLE*1e4*dt)*truthoffsetpartial{dt};
                        ydata = (xdata.*exp(-(xdata.^2)./(4*DMLE*frame_timedt))).*xdata+offset;
                        plot(xdata,(ydata./sum(ydata).*(size(xdata,2)/bgbinningnr)*Visualisation_HO_outputCell.bgratios(1)/Visualisation_HO_outputCell.bgratios(dt)),'k-','LineWidth',2)
                        %Get maxheight to later scale the subplots
                        maxheight(dt) = max(max(ydata),max(h.Values));
                    elseif populations == 2
                        %slightly different: plot both on top of baseline
                        %get baseline
                        yd2 = histInfo{dt}(:,2)*totalreconstit_arr_entries+baseline_to_be_added{dt};
                        baseline = ((baseline_to_be_added{dt}/mean(yd2))/size(histInfo{dt}(:,1),1))*ones(size(xdata));

                        DMLE1 = paramEsts(1)+loc_unc^2/frame_timedt;
                        offset1 = DMLE1*1e4*truthoffsetpartial{dt}*dt;
                        ydata1 = xdata.*exp(-(xdata.^2)./(4*DMLE1*frame_timedt)).*xdata+offset1;
                        ydata1norm = (ydata1./sum(ydata1).*(size(xdata,2)/bgbinningnr))*(1/(1+paramEsts(3)));
                        ydata1norm = ydata1norm-ydata1norm(end)+baseline;
                        plot(xdata,ydata1norm*Visualisation_HO_outputCell.bgratios(1)/Visualisation_HO_outputCell.bgratios(dt),'b-','LineWidth',1,'DisplayName',['MLE: ' num2str(paramEsts*1e12)]);

                        DMLE2 = paramEsts(2)+loc_unc^2/frame_timedt;
                        offset2 = DMLE2*1e4*truthoffsetpartial{dt}*dt;
                        ydata2 = xdata.*exp(-(xdata.^2)./(4*DMLE2*frame_timedt)).*xdata+offset2;
                        ydata2norm = (ydata2./sum(ydata2).*(size(xdata,2)/bgbinningnr))*(paramEsts(3)/(1+paramEsts(3)));
                        ydata2norm = ydata2norm-ydata2norm(end)+baseline;
                        plot(xdata,ydata2norm*Visualisation_HO_outputCell.bgratios(1)/Visualisation_HO_outputCell.bgratios(dt),'b-','LineWidth',1,'DisplayName',['MLE: ' num2str(paramEsts(1)*1e12)]);

                        ydatasum = ydata1norm-baseline+ydata2norm-baseline + baseline;
                        plot(xdata,ydatasum*Visualisation_HO_outputCell.bgratios(1)/Visualisation_HO_outputCell.bgratios(dt),'k-','LineWidth',1,'DisplayName',['MLE: ' num2str(paramEsts(2)*1e12) '-' num2str(paramEsts(3))]);

                        %Get maxheight to later scale the subplots
                        maxheight(dt) = max(max(ydatasum),max(h.Values));
                    end
                end
            else
                maxheight(dt) = max(h.Values);
            end
            %Plot the baseline
            yd2 = histInfo{dt}(:,2)*totalreconstit_arr_entries+baseline_to_be_added{dt};
            plot(histInfo{dt}(:,1),ones(size(histInfo{dt}(:,1)))./sum(ones(size(histInfo{dt}(:,1))))*(baseline_to_be_added{dt}/mean(yd2)),'k--')
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
        drawnow
    end
    Visualisation_HO_outputCell.baseline_to_be_added = baseline_to_be_added;
end % ------------ END OF APPROXIMATE FIT

%% Now start the 'original' MLE fit based on the data already obtained, with the variables obtained from HO
if performsecondfit
    if verbose; dispUIorCommandWindow('Starting full fit',callfromUI);end
    %Get inputarr list
    [size_dt,inputarr] = getSizesInputFromJDCell(JDarrSignalCell);
    pdfSettings.size_dt = size_dt;

    %BGcurve_interp gives a x,y list of interpolated points
    %     keyboard
    if strcmp(linorlogBGsubtract,'lin')
        BGcurve_interp = interpolate_BGCurve(maxdist,bgbinningnr,(maxdist/bgbinningnr),JDarrBG);
    else
        BGcurve_interp = interpolate_BGCurve_log(maxdist,minlogpoint,bgbinningnr,bgbinningnr,JDarrBG);
    end
    %output_BG_alldt is a list for all inputarr values on how likely it is that
    %it belongs to the BG. This is later rescaled according to a fraction of
    %the BG inside the MLE
    if strcmp(linorlogBGsubtract,'lin')
        output_BG_alldt = pdfBGFunction(inputarr,BGcurve_interp);
    else
        output_BG_alldt = pdfBGFunction_log(inputarr,BGcurve_interp);
    end

    %Set the MLEoptions
    mleoptions = TARDIS_setMLEoptions('FF');
    %%
    if settingsTARDIS.customPDF == 0
        if AutoChoosePop
            if verbose; dispUIorCommandWindow('Starting auto-determination of population nr...',callfromUI);end
            for nrpopsTesting = 1:2
                %Get the correct start, lower, upper parameters based on nr of
                %populations and estimated parameters etc
                pdfSettings.populations = nrpopsTesting;
                start = eval(eval(['start_' num2str(nrpopsTesting) 'pop']));
                eval(['lb = lb_' num2str(nrpopsTesting) 'pop;']);
                eval(['ub = ub_' num2str(nrpopsTesting) 'pop;']);
                clear startparam lowerbound upperbound
                [startparam, lowerbound, upperbound] = get_population_parametersFF(nrpopsTesting,paramEsts,HOfitCI,bgratios,start,lb_1pop,ub_1pop,lb_2pop,ub_2pop,performestimationfit);
    
                %Change parameters and lower/upper bound in case we fit with bleach
                %kinetics rather than fully random BG values
                [startparam, lowerbound, upperbound, fixRatios2pop] = changeFitParams_populations_Bleach(fitWithBleach, nrpopsTesting, startparam, lowerbound, upperbound, fixRatios2pop, dt_arr);
    
                %Create a customPDF
                %The used function can be found in the subscript defined here
                %(pdfBGwithPops), which has a certain set inputs, and the fitting variables
                %are in the varargin
    
                if strcmp(linorlogBGsubtract,'lin')
                    pdf_FF_pops = @(xdata,varargin)pdfBGwithPops(xdata,BGcurve_interp,output_BG_alldt,pdfSettings,callfromUI,varargin);
                else
                    pdf_FF_pops = @(xdata,varargin)pdfBGwithPops_log(xdata,BGcurve_interp,output_BG_alldt,pdfSettings,callfromUI,varargin);
                end
                %perform fit
                mleoptions.Display='off';
                [parameters_autoPickPops{nrpopsTesting}, parametersCI_autoPickPops{nrpopsTesting}] = mle(inputarr, 'pdf',pdf_FF_pops,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound,'Options',mleoptions);
    
                %Calculate adjusted R2 and store
                [R2_adj(nrpopsTesting), R2_adj_nopop(nrpopsTesting)] = R2_adj_calculation_TARDIS(maxdist, dt_arr, inputarr, bgbinningnr, JDarrBG, pdfSettings, callfromUI, parameters_autoPickPops{nrpopsTesting});
    
                if verbose; dispUIorCommandWindow(['Adjusted R-squared of ' num2str(nrpopsTesting) ' population(s): ' sprintf('%.6f', R2_adj(nrpopsTesting))],callfromUI);end
            end
            [~,bestNrPops] = max([R2_adj_nopop(1) R2_adj]);
            bestNrPops = bestNrPops-1; %In case 0-fit is best
            if verbose; dispUIorCommandWindow(['Best nr of populations: ' num2str(bestNrPops) ' - storing this info... '],callfromUI);end
            if bestNrPops == 0
                keyboard %To implement
            end
            %          keyboard
            parameters = parameters_autoPickPops{bestNrPops};
            parametersCI = parametersCI_autoPickPops{bestNrPops};
            populations = bestNrPops;
            pdfSettings.populations = bestNrPops;
        else %If auto-choose populations isn't used
            if populations > 0 %if single/multi population fit is wanted
                %Get the correct start, lower, upper parameters based on nr of
                %populations and estimated parameters etc
                [startparam, lowerbound, upperbound] = get_population_parametersFF(populations,paramEsts,HOfitCI,bgratios,start,lb_1pop,ub_1pop,lb_2pop,ub_2pop,performestimationfit);
    
                %Change parameters and lower/upper bound in case we fit with bleach
                %kinetics rather than fully random BG values
                [startparam, lowerbound, upperbound, fixRatios2pop] = changeFitParams_populations_Bleach(fitWithBleach, populations, startparam, lowerbound, upperbound, fixRatios2pop, dt_arr);
    
                %Create a customPDF
                %The used function can be found in the subscript defined here
                %(pdfBGwithPops), which has a certain set inputs, and the fitting variables
                %are in the varargin
                if strcmp(linorlogBGsubtract,'lin')
                    pdf_FF_pops = @(xdata,varargin)pdfBGwithPops(xdata,BGcurve_interp,output_BG_alldt,pdfSettings,callfromUI,varargin);
                else
                    pdf_FF_pops = @(xdata,varargin)pdfBGwithPops_log(xdata,BGcurve_interp,output_BG_alldt,pdfSettings,callfromUI,varargin);
                end
                %perform fit
                mleoptions.Display='off';
    
                %Change input for freefit loc unc
                % I leave this implemented, but MLE doesn't seem to work very
                % good
                % If freefit_locunc is here, add loc_unc to the fitting params
                if pdfSettings.freefit_locunc
                    startparam = [startparam loc_unc*1e9];
                    lowerbound = [lowerbound max(loc_unc*1e9-50,5)];
                    upperbound = [upperbound loc_unc*1e9+50];
                end
                %Fix small error if startparam is ever same value as
                %upperbound/lowerbound
                for i = 1:size(startparam,2)
                    if startparam(i)<=lowerbound(i)
                        startparam(i) = lowerbound(i)*1.0001;
                        disp('CHANGED MLE FIT STARTING PARAMS PROBABLY SOMETHING WRONG')
                    end
                    if startparam(i)>=upperbound(i)
                        startparam(i) = upperbound(i)*0.9999;
                        disp('CHANGED MLE FIT STARTING PARAMS PROBABLY SOMETHING WRONG')
                    end
                end
                [parameters, parametersCI] = mle(inputarr, 'pdf',pdf_FF_pops,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound,'Options',mleoptions);
                %After fitting, extract the loc_unc again and restore the
                %parameters to what's expected. Also report on loc_unc
                %confidence interval, which is normally pretty bad.
                if pdfSettings.freefit_locunc
                    [parameters, parametersCI] = mle(inputarr, 'pdf',pdf_FF_pops,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound,'Options',mleoptions);
                    loc_unc = parameters(end)/1e9;
                    disp(['Found Loc unc: ' num2str(loc_unc*1e9) '(' num2str(parametersCI(1,end)) ' - ' num2str(parametersCI(2,end)) ')'])
                    parameters = parameters(1:end-1);
                    parametersCI = parametersCI(:,1:end-1);
                end
                %%
                %Calculate R^2
                %             [R2_adj, R2_adj_nopop] = R2_adj_calculation_TARDIS(maxdist, dt_arr, inputarr, bgbinningnr, JDarrBG, frame_time, loc_unc, size_dt, populations, fitWithBleach, fixRatios2pop, callfromUI, parameters);
                [R2_adj, R2_adj_nopop] = R2_adj_calculation_TARDIS(maxdist, dt_arr, inputarr, bgbinningnr, JDarrBG, pdfSettings, callfromUI, parameters);
                %%
                if populations == 1
                    if verbose
                        stringval = sprintf('Second fit (final) completed with DiffCoeff %.2f (%.2f-%.2f).\n',parameters(1),parametersCI(1,1),parametersCI(2,1));
                        dispUIorCommandWindow(stringval,callfromUI);
                    end
                elseif populations == 2
                    if verbose
                        stringval = sprintf('Second fit (final) completed with DiffCoeffs: %.2f (%.2f-%.2f) and %.2f (%.2f-%.2f) at ratio 1:%.2f (%.2f-%.2f).\n',parameters(1),parametersCI(1,1),parametersCI(2,1),parameters(2),parametersCI(1,2),parametersCI(2,2),parameters(3),parametersCI(1,3),parametersCI(2,3));
                        dispUIorCommandWindow(stringval,callfromUI);
                    end
                end
            end
        end
        if AutoChoosePop == 0 && populations == 0 %anaDDA fit - never auto-choose nr of pops
            %Create anaDDA input - following anaDDA matlab
            % Fitting function
    %         keyboard
            %Obtain good BG_alldt_D list (with BG information)
            for dt = 1:size(size_dt,2)
                output_BG_alldt_D{dt} = pdfBGFunction_log(inputarr_D,BGcurve_interp_D{dt});
            end
    
            % %Set the MLEoptions
            mleoptions = TARDIS_setMLEoptions('FF');
    
            if performestimationfit == 0
                startparam = [str2num(settingsTARDIS.start_aDDA) fliplr(bgratios)]';
                paramEsts = [0 0 0];
                HOfitCI = [0 0 0];
            else
                startparam = [paramEsts(1:3)' fliplr(bgratios)]';
            end
            lowerbound = [lb_aDDA ones(1,max(dt_arr))*0]';
            upperbound = [ub_aDDA ones(1,max(dt_arr))*1]';
    
            %Create anaDDA custom PDF
            pdf_FF_aDDA = @(Dlistdata,varargin) pdfAnaDDAMLE_multidt(inputarr_D,BGcurve_interp_D, Numberofframes,Frametimelist,...
                inputAnaDDA,rangeD_alldtbins,fitspecies_singledtbin{1}, fixedspecies_singledtbin{1}, ...
                (((Dfixed_singledtbin{dt}(1)))), fixedparameters,indexfittingparameters,fx,fy,maxDindtracking,...
                frequency,size_dt,(locerrorpdfcorrected_singledtbin),output_BG_alldt_D,bgbinningnr,fitWithBleach,verbose,callfromUI,varargin);
    
    %         keyboard
    %             %%
    %             figure(92);clf(92);
    %             histogram(inputarr_D(1:size_dt(1)),logspace((-4),log10(3),100))
    %             hold on
    %             plot(BGcurve_interp_D{1})
    %             set(gca,'XScale','log')
            %%
            %Change parameters and lower/upper bound in case we fit with bleach
            %kinetics rather than fully random BG values
            if fitWithBleach
                [startparam, lowerbound, upperbound] = changeFitParamsBleach([0.08,0.0001,0.5],startparam,lowerbound,upperbound,'aDDA');
            end
    
            %Perform anaDDA MLE fit
            [parameters] = mle(Dlistdata, 'pdf',pdf_FF_aDDA,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound,'Options',mleoptions);
            %Calculate convindence interval via coveriance matrix
            acov = mlecov(parameters,Dlistdata,'pdf',pdf_FF_aDDA);
            se = sqrt(diag(acov))';
            %Calculate 95pct confidence interval
            alpha = 1-0.95;
            probs = [alpha/2; 1-alpha/2];
            phat = parameters';
            pci = norminv(repmat(probs,1,numel(phat)),[phat; phat],[se; se]);
            parametersCI = pci;
        end
    else
        %% Custom PDF
        [parameters, parametersCI] = TARDIS_customPDF_wrapper(customPDF, inputarr, bgratios, linorlogBGsubtract, BGcurve_interp,output_BG_alldt,pdfSettings, eval(customPDF_start), eval(customPDF_lower), eval(customPDF_upper), callfromUI);
    end
    %% Visualise results from MLE
    if populations > 0
        figurebinnr = settingsTARDIS.bgbinningnr;
        maxdtbins = max(dt_arr);
        xdataplotsingledt = [linspace(0,maxdist,(figurebinnr+1))]'; %Create bins used for visualisation, based on figurebinnr
        xdataplot = repmat(xdataplotsingledt,maxdtbins,1); %Extend these bins for every dtbin, and put them below each other
        size_dt = ones(1,maxdtbins)*(figurebinnr+1); %Make a 1-by-dt array, each holding the size of the bins, to later split up xdataplot again

        %Always linear BGcurve interp for visualisation
        BGcurve_interp_vis = interpolate_BGCurve(maxdist,bgbinningnr,(maxdist/bgbinningnr),JDarrBG);
        output_BG_alldt = pdfBGFunction(xdataplot,BGcurve_interp_vis); %Get BG values for all x positions (repeating for dt positions)
        %Get signal values for all x positions and all dt positions
        pdfSettings.verbose = 0;
        pdfSettings.size_dt = size_dt;
        [~,extraoutput] = pdfBGwithPops(xdataplot,BGcurve_interp_vis,output_BG_alldt,pdfSettings,callfromUI,{parameters});
        %Save visualisation output for later out-script use
        visualisationinfo.size_dt = size_dt;
        visualisationinfo.JDarrSignalCell = JDarrSignalCell;
        visualisationinfo.xdataplot = xdataplot;
        visualisationinfo.BGcurve_interp = BGcurve_interp_vis;
        visualisationinfo.extraoutput = extraoutput;
        visualisationinfo.frame_time = frame_time;
        visualisationinfo.fixRatios2pop = fixRatios2pop;
        % %         %Start actual visualistion
        % %         if visualisationMLEIntFit
        % %             Visualisation_FF('','',size_dt,JDarrSignalCell,xdataplot,BGcurve_interp_vis,parameters,fitWithBleach,populations,'lin',visualisationinfo);
        % %         end

        Visualisation_FF_outputCell.size_dt = size_dt;
        Visualisation_FF_outputCell.poslist = poslist;
        Visualisation_FF_outputCell.JDarrSignalCell = JDarrSignalCell;
        Visualisation_FF_outputCell.xdataplot = xdataplot;
        Visualisation_FF_outputCell.BGcurve_interp = BGcurve_interp_vis;
        Visualisation_FF_outputCell.FFparameters = parameters;
        Visualisation_FF_outputCell.populations = populations;
        Visualisation_FF_outputCell.fitWithBleach = fitWithBleach;
        Visualisation_FF_outputCell.fixRatios2pop = fixRatios2pop;
        Visualisation_FF_outputCell.linorlog = settingsTARDIS.linorlogvis;
        Visualisation_FF_outputCell.strobo_frame_time = pdfSettings.strobo_frame_time;
        Visualisation_FF_outputCell.extraoutput.main = extraoutput;
        Visualisation_FF_outputCell.extraoutput.loc_unc = loc_unc;
        Visualisation_FF_outputCell.extraoutput.frame_time = frame_time;
        Visualisation_FF_outputCell.extraoutput.strobo_frame_time = pdfSettings.strobo_frame_time;
        Visualisation_FF_outputCell.extraoutput.fixRatios2pop = fixRatios2pop;
        Visualisation_FF_outputCell.extraoutput.logvis.minlogval = minlogpoint;
        Visualisation_FF_outputCell.extraoutput.logvis.BGcurve_interplog = interpolate_BGCurve_log(maxdist,minlogpoint,bgbinningnr,bgbinningnr,JDarrBG);
        Visualisation_FF_outputCell.extraoutput.logvis.xdataplotlog = repmat(logspace(log10(minlogpoint),log10(maxdist),bgbinningnr+1),1,size(size_dt,2))';

        %Start actual visualistion
        if visualisationMLEIntFit
            Visualisation_FF('','',Visualisation_FF_outputCell);
            %             Visualisation_FF('','',size_dt,JDarrSignalCell,xdataplot,BGcurve_interp_vis,parameters,fitWithBleach,populations,'log',Visualisation_FF_outputCell.extraoutput);
        end

    elseif customPDF == 0 %anaDDA
        [~,anaDDAvisInfoFF.visual] = pdfAnaDDAMLE_multidt(inputarr_D,BGcurve_interp_D, Numberofframes,ones(size(inputarr_D))*frame_time,...
            inputAnaDDA,rangeD_alldtbins,fitspecies_singledtbin{1}, fixedspecies_singledtbin{1}, ...
            (((Dfixed_singledtbin{dt}(1)))), fixedparameters,indexfittingparameters,fx,fy,maxDindtracking,...
            frequency,size_dt,(locerrorpdfcorrected_singledtbin),output_BG_alldt_D,bgbinningnr,fitWithBleach,verbose,callfromUI,parameters);

        %Honeslty no idea why I change settings here:
        figurebinnr = settingsTARDIS.bgbinningnr;
        maxdtbins = max(dt_arr);
        xdataplotsingledt = [linspace(0,maxdist,(figurebinnr+1))]'; %Create bins used for visualisation, based on figurebinnr
        xdataplot = repmat(xdataplotsingledt,maxdtbins,1); %Extend these bins for every dtbin, and put them below each other
        size_dt = ones(1,maxdtbins)*(figurebinnr+1); %Make a 1-by-dt array, each holding the size of the bins, to later split up xdataplot again
        %output_BG_alldt = pdfBGFunction_log(xdataplot,BGcurve_interp); %Get BG values for all x positions (repeating for dt positions)

        %         output_BG_alldt = pdfBGFunction(xdataplot,BGcurve_interp); %Get BG values for all x positions (repeating for dt positions)
        %         %Get signal values for all x positions and all dt positions
        %
        % %         keyboard
        %         %Visualisation info for HO-anaDDA
        % %                 keyboard
        % % keyboard
        %         [~,anaDDAvisInfoFF.visual] = pdfAnaDDAMLE_multidt(inputReconstitarr_D,BGcurve_interp_D, Numberofframes_reconstit,Frametimelist,...
        %             inputAnaDDA,rangeD_alldtbins,fitspecies_singledtbin{1}, fixedspecies_singledtbin{1}, ...
        %             (((Dfixed_singledtbin{dt}(1)))), fixedparameters,indexfittingparameters,fx,fy,maxDindtracking,...
        %             frequency_reconstit,size_dt_reconstit,(locerrorpdfcorrected_singledtbin),output_BG_alldt_D_linear,bgbinningnr,fitWithBleach,verbose,callfromUI,parameters);

        anaDDAvisInfoFF.DlistFF = inputarr_D;
        anaDDAvisInfoFF.DlistHO = inputReconstitarr_D;
        anaDDAvisInfoFF.size_dt_FF = frequency;
        anaDDAvisInfoFF.size_dt_HO = frequency_reconstit;
        %Save visualisation output for later out-script use
        visualisationinfo.size_dt = size_dt;
        visualisationinfo.JDarrSignalCell = JDarrSignalCell;
        visualisationinfo.xdataplot = xdataplot;
        visualisationinfo.BGcurve_interp = BGcurve_interp_D;
        anaDDAvisInfoFF.visualisationinfo = visualisationinfo;
        anaDDAvisInfoFF.main = anaDDAvisInfoFF;
        visualisationinfo.extraoutput = anaDDAvisInfoFF;
        visualisationinfo.fitWithBleach = fitWithBleach;

        Visualisation_FF_outputCell.size_dt = size_dt;
        Visualisation_FF_outputCell.JDarrSignalCell = JDarrSignalCell;
        Visualisation_FF_outputCell.xdataplot = xdataplot;
        Visualisation_FF_outputCell.BGcurve_interp = BGcurve_interp;
        Visualisation_FF_outputCell.FFparameters = parameters;
        Visualisation_FF_outputCell.fitWithBleach = fitWithBleach;
        Visualisation_FF_outputCell.populations = populations;
        Visualisation_FF_outputCell.fixRatios2pop = fixRatios2pop;
        Visualisation_FF_outputCell.linorlog = settingsTARDIS.linorlogvis;
        %         Visualisation_FF_outputCell.extraoutput = anaDDAvisInfoFF;
        Visualisation_FF_outputCell.extraoutput.main = anaDDAvisInfoFF;

        Visualisation_FF_outputCell.extraoutput.loc_unc = loc_unc;
        Visualisation_FF_outputCell.extraoutput.frame_time = frame_time;
        Visualisation_FF_outputCell.extraoutput.strobo_frame_time = pdfSettings.strobo_frame_time;
        Visualisation_FF_outputCell.extraoutput.logvis.minlogval = minlogpoint;
        Visualisation_FF_outputCell.extraoutput.linvis.maxdist = maxdist;
        Visualisation_FF_outputCell.extraoutput.linvis.BGarr_maxdist_D = BGarr_maxdist_D;
        for dt = 1:size(size_dt,2)
            Visualisation_FF_outputCell.extraoutput.linvis.BGcurve_interp_lin{dt} = interpolate_BGCurve(maxdist_D{dt},bgbinningnr,maxdist_D{dt}./bgbinningnr,BGarr_maxdist_D{dt});
        end

        %Start actual visualistion
        if visualisationMLEIntFit
            %             Visualisation_FF('','',size_dt,JDarrSignalCell,xdataplot,BGcurve_interp,parameters,fitWithBleach,populations,anaDDAvisInfo);
            Visualisation_FF('','',Visualisation_FF_outputCell);
            %             Visualisation_FF('','',size_dt,JDarrSignalCell,xdataplot,BGcurve_interp,parameters,fitWithBleach,populations,linorlogvis,anaDDAvisInfoFF);
        end
        %         Visualisation_FF_outputCell.extraoutput.linvis.BGcurve_interplin_D = interpolate_BGCurve(maxdist,minlogpoint,bgbinningnr,bgbinningnr,JDarrBG);
        %         Visualisation_FF_outputCell.extraoutput.logvis.xdataplotlin = repmat(logspace(log10(minlogpoint),log10(maxdist),bgbinningnr+1),1,size(size_dt,2))';

    elseif customPDF ~= 0
        %Need a few empty parameters as output
        paramEsts=[];
        HOfitCI = [];

        %Probably need to change something here for visualisation down the
        %line
        figurebinnr = settingsTARDIS.bgbinningnr;
        maxdtbins = max(dt_arr);
        xdataplotsingledt = [linspace(0,maxdist,(figurebinnr+1))]'; %Create bins used for visualisation, based on figurebinnr
        xdataplot = repmat(xdataplotsingledt,maxdtbins,1); %Extend these bins for every dtbin, and put them below each other
        size_dt = ones(1,maxdtbins)*(figurebinnr+1); %Make a 1-by-dt array, each holding the size of the bins, to later split up xdataplot again

        %Always linear BGcurve interp for visualisation
        BGcurve_interp_vis = interpolate_BGCurve(maxdist,bgbinningnr,(maxdist/bgbinningnr),JDarrBG);
        output_BG_alldt = pdfBGFunction(xdataplot,BGcurve_interp_vis); %Get BG values for all x positions (repeating for dt positions)
        %Get signal values for all x positions and all dt positions
        pdfSettings.verbose = 0;
        pdfSettings.size_dt = size_dt;
        [~,extraoutput] =  CustomPDFwithBG(xdataplot,BGcurve_interp_vis,output_BG_alldt,pdfSettings,callfromUI,customPDF,{parameters});
%         [~,extraoutput] = pdfBGwithPops(xdataplot,BGcurve_interp_vis,output_BG_alldt,pdfSettings,callfromUI,{parameters});
        %Save visualisation output for later out-script use
        visualisationinfo.size_dt = size_dt;
        visualisationinfo.JDarrSignalCell = JDarrSignalCell;
        visualisationinfo.xdataplot = xdataplot;
        visualisationinfo.BGcurve_interp = BGcurve_interp_vis;
        visualisationinfo.extraoutput = extraoutput;
        visualisationinfo.frame_time = frame_time;
        visualisationinfo.fixRatios2pop = fixRatios2pop;

        Visualisation_FF_outputCell.size_dt = size_dt;
        Visualisation_FF_outputCell.poslist = poslist;
        Visualisation_FF_outputCell.JDarrSignalCell = JDarrSignalCell;
        Visualisation_FF_outputCell.xdataplot = xdataplot;
        Visualisation_FF_outputCell.BGcurve_interp = BGcurve_interp_vis;
        Visualisation_FF_outputCell.FFparameters = parameters;
        Visualisation_FF_outputCell.populations = populations;
        Visualisation_FF_outputCell.fitWithBleach = fitWithBleach;
        Visualisation_FF_outputCell.fixRatios2pop = fixRatios2pop;
        Visualisation_FF_outputCell.linorlog = settingsTARDIS.linorlogvis;
        Visualisation_FF_outputCell.strobo_frame_time = pdfSettings.strobo_frame_time;
        Visualisation_FF_outputCell.extraoutput.main = extraoutput;
        Visualisation_FF_outputCell.extraoutput.loc_unc = loc_unc;
        Visualisation_FF_outputCell.extraoutput.frame_time = frame_time;
        Visualisation_FF_outputCell.extraoutput.strobo_frame_time = pdfSettings.strobo_frame_time;
        Visualisation_FF_outputCell.extraoutput.fixRatios2pop = fixRatios2pop;
        Visualisation_FF_outputCell.extraoutput.logvis.minlogval = minlogpoint;
        Visualisation_FF_outputCell.extraoutput.logvis.BGcurve_interplog = interpolate_BGCurve_log(maxdist,minlogpoint,bgbinningnr,bgbinningnr,JDarrBG);
        Visualisation_FF_outputCell.extraoutput.logvis.xdataplotlog = repmat(logspace(log10(minlogpoint),log10(maxdist),bgbinningnr+1),1,size(size_dt,2))';

        Visualisation_FF_outputCell.customPDF = customPDF;

        %Start actual visualistion
        if visualisationMLEIntFit
            Visualisation_FF_customPDF('','',Visualisation_FF_outputCell);
            %             Visualisation_FF('','',size_dt,JDarrSignalCell,xdataplot,BGcurve_interp_vis,parameters,fitWithBleach,populations,'log',Visualisation_FF_outputCell.extraoutput);
        end
%         Visualisation_FF_outputCell=[];
    end
else
    disp('No second fit being performed (settings.performsecondfit set to false)')
    parameters=[];
end
%% Finalize
tottime = toc;

dispUIorCommandWindow('Finished fitting!',callfromUI);
dispUIorCommandWindow('-----------------------------------------',callfromUI);
%Two more enters only when from command line
if isempty(callfromUI)
    dispUIorCommandWindow(' ',callfromUI);
    dispUIorCommandWindow(' ',callfromUI);
end
end