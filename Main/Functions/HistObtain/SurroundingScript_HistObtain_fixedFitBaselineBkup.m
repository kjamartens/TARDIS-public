load('Z:\Documents\WUR info - data\Backup H Drive\MATLAB_BIPNAS\RelativeDistance_trackingAlternative\SimulatedTracks\OnePop_anaDDA\dens_2_kon_0_koff_0_Dfree_100.mat');
tic
frame_dist_BG = 10;
maxdist=2e-6;
frame_dist_Hist_arr=[1];
vis = 0;
bgbinningnr = 100;
[output,time,bgarr] = MLE_BG_subtraction_HistObtain_function(poslist,frame_dist_BG,maxdist,frame_dist_Hist_arr,bgbinningnr,vis,poslistunordered);



% LSQFIT
frame_time = 0.01;
loc_unc = 15e-9*1e6;
dtbins = frame_dist_Hist_arr;
fun = @(D,xdata)LSQfitter(xdata,D,loc_unc,frame_time,dtbins);

startD = [2e-12]*1e12; 
clear options
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','OptimalityTolerance',1e-12,'FunctionTolerance',1e-9);%,'StepTolerance',1e-15);
% options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','OptimalityTolerance',1e-12);%,'FunctionTolerance',1e-9);%,'StepTolerance',1e-15);
lb = [0.01e-12]*1e12;
ub = [10e-12]*1e12;

fitxdata = repmat(output{1}(:,1)*1e6,size(dtbins,2),1);
fitydata = zeros(size(fitxdata));
for dt = 1:size(dtbins,2)
    startp = 1+size(output{1}(:,2),1)*(dt-1);
    endp = startp+size(output{1}(:,2),1)-1;
    fitydata(startp:endp) = output{dt}(:,2);
end
clear startp endp
outp = lsqcurvefit(fun,startD,fitxdata,fitydata,lb,ub,options);
toc
outp

%% 
figure(4);clf(4);
for dt = 1:3
subplot(1,3,dt);
hold on
plot(output{1}(:,1),output{dt}(:,2),'r-x')
xdat = [0:1e-8:2e-6];
% Df = outp/1e12;
Df = 1/1e12;
loc_unc = 15e-9;
frt = frame_time*dt;
ydat = xdat.*exp(-(xdat.^2)./(4*(Df+loc_unc^2/frt)*frt));
ydat = ydat./sum(ydat)*2;
plot(xdat,ydat,'k-')
axis([0 1e-6 0 0.1])
end
%% Create 'fake' JD array
%it might have slightly decreased or biased results because of removal of
%negatives

%find edges
halfwidthsize = (output{1}(2,1)-output{1}(1,1))/2;
if (halfwidthsize-output{1}(1,1)) > 1e-16
    disp('This is going to fail! Probably change code')
end
startpoint = 0;
endpoint = maxdist;
edges = [startpoint:halfwidthsize*2:endpoint];

totalreconstit_arr_entries = 10000; % This is a sort of minimum? It's slightly higher than this because of removal of negatives
baseline_per_bin = true;
reconstit_arr = zeros(totalreconstit_arr_entries,1);
counter = 1;
summer = [];
baseline_to_be_added = 0;
if baseline_per_bin
    for i = 1:size(edges,2)-1
    %find nr of entries per bin
    nrentries_bin_check(i) = floor(output{1}(i,2)*totalreconstit_arr_entries);
    
    end
    %get minimum and add later
    baseline_to_be_added = max(0,min(nrentries_bin_check)*-1+1);
end
for i = 1:size(edges,2)-1
    if ~baseline_per_bin
        nrentries_bin = max(round(output{1}(i,2)*totalreconstit_arr_entries),0);
    else
        nrentries_bin = round(output{1}(i,2)*totalreconstit_arr_entries)+baseline_to_be_added;
    end
    summer(i) = nrentries_bin;
    reconstit_arr(counter:counter+nrentries_bin-1) = rand(nrentries_bin,1)*(halfwidthsize*2)+edges(i);
%     reconstit_arr(counter:counter+nrentries_bin-1) = (halfwidthsize)+edges(i);

    counter = counter + nrentries_bin;
end

%Next value is the ratio of numbers that are there for the baseline rather
%than for the data.
ratio_baseline_to_data = (baseline_to_be_added*(size(edges,2)-1))/size(reconstit_arr,1);


% MLE fit again

frame_time = 0.01;
loc_unc = 15e-9;
norm_bins = 10000;

truthoffset = 1e-8*(baseline_to_be_added*bgbinningnr)/(size(reconstit_arr,1)-(baseline_to_be_added*bgbinningnr))

custompdf = @(x,D)MLEfitter(x,D,truthoffset,loc_unc,frame_time,dtbins,norm_bins,ratio_baseline_to_data);
% custompdf = @(x,D)MLEfitter(x,D,loc_unc,frame_time,dtbins,norm_bins,ratio_baseline_to_data);
reconstit_arr_MLE = reconstit_arr;%(reconstit_arr(reconstit_arr<1e-6));
start = [20e-12];
lb = [1e-15];
ub = [30e-12];
% start = [1e-12];
% lb = [1e-15];
% ub = [30e-12];
mleoptions = statset('mlecustom');
mleoptions.MaxIter = 2000; %Maximum iterations
mleoptions.MaxFunEvals = mleoptions.MaxIter; %Maximum function evaluations
tic
paramEsts = mle(reconstit_arr_MLE, 'pdf',custompdf, 'start',start, ...
    'lower',lb, 'upper',ub, 'options',mleoptions)
paramEsts(1)*1e12
% paramEsts(2)*1e12
toc
%% Visualisation
figure(3);clf(3);
for dt = frame_dist_Hist_arr
subplot(1,3*max(frame_dist_Hist_arr),1+3*(dt-1))
histogram(reconstit_arr,linspace(0,maxdist,bgbinningnr+1));%,'Normalization','Probability');
hold on
xdata = [0:1e-8:2e-6];
plot(xdata,ones(size(xdata))*baseline_to_be_added,'r-');
origdataoutput = output{1}(:,2).*(size(reconstit_arr,1)-baseline_to_be_added*(size(edges,2)-1));
origdataoutput = origdataoutput+baseline_to_be_added;
plot(output{1}(:,1),origdataoutput,'m-x')

subplot(1,3*max(frame_dist_Hist_arr),2+3*(dt-1))
origdataoutput = output{1}(:,2);
plot(output{1}(:,1),origdataoutput,'m-x')
hold on

loc_unc = 15e-9;
frame_time = 0.01;
DMLEt = 1e-12+loc_unc^2/frame_time;
ydata = output{1}(:,1).*exp(-(output{1}(:,1).^2)./(4*DMLEt*frame_time));
plot(output{1}(:,1),ydata./sum(ydata),'b-','LineWidth',2,'DisplayName',['MLE: ' num2str(paramEsts*1e12)]);


subplot(1,3*max(frame_dist_Hist_arr),3+3*(dt-1))
histogram(reconstit_arr,100,'Normalization','Probability');
hold on
xdata = [0:1e-8:2e-6];
loc_unc = 15e-9;
frame_time = 0.01;
% flatline = (ratio_baseline_to_data.*size(reconstit_arr,1))./1e12;
flatline = ratio_baseline_to_data;
MLEt = 1e-12+loc_unc^2/frame_time;
yd1 = xdata.*exp(-(xdata.^2)./(4*MLEt*frame_time));
yd1 = yd1./sum(yd1);
yd2 = ones(size(yd1));
yd2 = yd2./sum(yd2);
ydata = yd1.*(1-ratio_baseline_to_data)+yd2.*(ratio_baseline_to_data);
ydata = ydata./sum(ydata);
plot(xdata,ydata.*(size(xdata,2)/100),'k-','LineWidth',2,'DisplayName',['MLE: ' num2str(paramEsts*1e12)]);

DMLE = paramEsts(1);
DMLE=DMLE+loc_unc^2/frame_time;
ydata = xdata.*exp(-(xdata.^2)./(4*DMLE*frame_time))+truthoffset;
% ydd1 = xdata.*exp(-(xdata.^2)./(4*DMLE*frame_time));
% ydd2 = truthoffset;
% sum(ydd2)/sum(ydd1)
% sum(ydd2./sum(ydata).*(size(xdata,2)/100))./sum(ydd1./sum(ydata).*(size(xdata,2)/100))
plot(xdata,ydata./sum(ydata).*(size(xdata,2)/100),'c--','LineWidth',2,'DisplayName',['MLE: ' num2str(paramEsts*1e12)]);


DLSQ = outp/1e12+loc_unc^2/frame_time;
% ydat = xdata.*exp(-(xdata.^2)./(4*(DLSQ+loc_unc^2/frt)*frt));
% ydat = ydat./sum(ydat);
% plot(xdata,ydat./sum(ydat),'b-','LineWidth',2,'DisplayName',['LSQ: ' num2str(outp)]);

axis([0 2e-6 0 inf])
% legend()

end
%% Poor-mans fitting

%
clear output_Diff
Darr = [0.1e-12:0.001e-12:2e-12];
for Dcntr = 1:size(Darr,2)
    D = Darr(Dcntr);
    loc_unc = 15e-9;
    frame_time = 0.01;
    D=D+loc_unc^2/frame_time;
    xdata = output{1}(:,1);
    output_Diff(Dcntr,:) =  xdata.*exp(-(xdata.^2)./(4*D*frame_time));%exp(-(xdata./(4*D*frame_time)));
    output_Diff(Dcntr,:) = output_Diff(Dcntr,:)./sum(output_Diff(Dcntr,:));
    rmse(Dcntr) = sum(sqrt((output_Diff(Dcntr,:)'-output{1}(:,2)).^2));
end
[a,b] = min(rmse);
minDfound = Darr(b)
% output_Diff = output_Diff./(size(output{1},2)/size(output_Diff,2));



%%
figure(2);clf(2);
subplot(1,2,1)
plot(output{1}(:,1),output{1}(:,2)*10000/size(reconstit_arr,1),'k-x')
axis([0 1e-6 0 inf])
hold on
histogram(reconstit_arr,linspace(0,maxdist,bgbinningnr+1),'Normalization','Probability');
% plot(xdata,output_Diff(find(abs(Darr-1e-12)<1e-16),:),'b-') %correct one
% plot(xdata,output_Diff(b,:),'r-') % found one
subplot(1,2,2)
plot(bgarr(:,1),bgarr(:,2),'k-x')

%% MLE fit of multiple dt
frame_time = 0.01;
populations = 1;
size_dt = [11059];
verbose = 1;
BGcurve_interp = bgarr;
inputarr = reconstit_arr;


%% Testing two univeriates
% clc;
% custompdf = @(x,varargin)MLEfitter(x,frame_time,populations,size_dt,verbose,varargin);
% inputarr = (reconstit_arr(reconstit_arr<0.5e-6));
% start = [0.8e-12];
% lb = [1e-15];
% ub = [2e-12];
% options = statset('MaxIter',500, 'MaxFunEvals',6000);
% paramEsts = mle(inputarr, 'pdf',custompdf, 'start',start, ...
%     'lower',lb, 'upper',ub, 'options',options)
%%
%%

%%
% function [output] = MLEfitter(inputdata,frame_time,populations,size_dt,verbose,varargin)
% vars = cell2mat(varargin{1,1});
% D = vars(1)*1e12;
% loc_unc = 15e-9;
% frame_time = 0.01;
% D=D+loc_unc^2/frame_time;
% % output =  inputdata.*exp(-(inputdata.^2)./(4*D*frame_time));
% % output =  exp(-(inputdata.^2)./(4*D*frame_time));
% output = normpdf(inputdata,D,D/10);
% end

%% function LSQ
function [output] = LSQfitter(xdata,D,loc_unc,frame_time,dtbins)
%assuming that f.e. dtbins = [1 2 3], meaning 1,2,3 frame delay, and all
%same x-range/xdata
%Output should be [ydatadt1 ydatadt2 ydatadt3]
% keyboard
xdatasingle = xdata(1:size(xdata,1)/size(dtbins,2));

for dt = dtbins
    psuedoframetime = dt.*frame_time;
    outputsingledt{dt} = xdatasingle.*exp(-(xdatasingle.^2)./(4*(D+loc_unc^2/psuedoframetime)*psuedoframetime));
    outputsingledt{dt} = outputsingledt{dt}./sum(outputsingledt{dt});
end

output = zeros(size(xdata));
for dt = 1:size(dtbins,2)
    startp = 1+size(outputsingledt{1},1)*(dt-1)
    endp = startp+size(outputsingledt{1},1)-1
    output(startp:endp) = outputsingledt{dt};
end

end

%% function MLE
function [output] = MLEfitter(xdata,D,offset,loc_unc,frame_time,dtbins,norm_bins,ratio_baseline_to_data)
% function [output] = MLEfitter(xdata,D,loc_unc,frame_time,dtbins,norm_bins,ratio_baseline_to_data)
%assuming that f.e. dtbins = [1 2 3], meaning 1,2,3 frame delay, and all
%same x-range/xdata
%Output should be [ydatadt1 ydatadt2 ydatadt3]
% keyboard
% keyboard
% output = pdfSinglePopFunction(xdata,D,offset,loc_unc,frame_time);
% output = normalize_area(xdata,output,D,offset,loc_unc,frame_time,norm_bins);
% offset = 0;

%This is sorta working
output = pdfSinglePopFunction(xdata,D,offset,loc_unc,frame_time);
output = normalize_area(xdata,output,D,offset,loc_unc,frame_time,norm_bins);
output(output<=0) = 1e-150;

% ratio_baseline_to_data
%this is experimental and should be worked on
% offset = 0;
% output = pdfSinglePopFunction_offsetbaseline(xdata,D,offset,loc_unc,frame_time,ratio_baseline_to_data);
% output = normalize_area_offsetbaseline(xdata,output,D,offset,loc_unc,frame_time,norm_bins,ratio_baseline_to_data);
% output(output<=0) = 1e-150;
end
%%
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
%%

function [output] = normalize_area_offsetbaseline(xdata,input,D,offset,loc_unc,frame_time,norm_bins,ratio_baseline_to_data)
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
mean_y_pos_bin = pdfSinglePopFunction_offsetbaseline(mean_x_pos_bin,D,offset,loc_unc,frame_time,ratio_baseline_to_data);
%Get total area
area_bin = mean_y_pos_bin*area_calc_bin_size;
totarea = sum(area_bin);
%Correct for total area
output = input./totarea;
end

%%
function output = pdfSinglePopFunction(xdata,D,offset,loc_unc,frame_time)
    output = xdata.*exp(-(xdata.^2)./(4*(D+loc_unc^2/frame_time)*frame_time))+offset;
end
%%
function output = pdfSinglePopFunction_offsetbaseline(xdata,D,offset,loc_unc,frame_time,ratio_baseline_to_data)
%     output = xdata.*exp(-(xdata.^2)./(4*(D+loc_unc^2/frame_time)*frame_time))+offset;

    yd1 = xdata.*exp(-(xdata.^2)./(4*(D+loc_unc^2/frame_time)*frame_time));
    yd1 = yd1./sum(yd1);
    yd2 = ones(size(yd1));
    yd2 = yd2./sum(yd2);
    output = yd1.*(1-ratio_baseline_to_data)+yd2.*(ratio_baseline_to_data);
end