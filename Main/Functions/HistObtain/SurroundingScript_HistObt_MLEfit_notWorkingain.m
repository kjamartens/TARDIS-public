load('Z:\Documents\WUR info - data\Backup H Drive\MATLAB_BIPNAS\RelativeDistance_trackingAlternative\SimulatedTracks\OnePop_anaDDA\dens_3_kon_0_koff_0_Dfree_100.mat');
frame_dist_BG = 10;
maxdist=2e-6;
frame_dist_Hist_arr=[1 2 3];
vis = 0;
bgbinningnr = 100;
[output,time,bgarr] = MLE_BG_subtraction_HistObtain_function(poslist,frame_dist_BG,maxdist,frame_dist_Hist_arr,bgbinningnr,vis,poslistunordered);
time

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
reconstit_arr = zeros(totalreconstit_arr_entries,1);
counter = 1;
summer = [];
for i = 1:size(edges,2)-1
    nrentries_bin = max(round(output{1}(i,2)*totalreconstit_arr_entries),0);
    summer(i) = nrentries_bin;
    reconstit_arr(counter:counter+nrentries_bin-1) = rand(nrentries_bin,1)*(halfwidthsize*2)+edges(i);
%     reconstit_arr(counter:counter+nrentries_bin-1) = (halfwidthsize)+edges(i);
    
    counter = counter + nrentries_bin;
end

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
%%
% output_BG_alldt = pdfBGFunction(inputarr,BGcurve_interp);
%This is done on a (fake) JD array for future integration with anaDDA
custompdf = @(xdata,varargin)MLEfitter(xdata,frame_time,populations,size_dt,verbose,varargin);
% custompdf = @(xdata,varargin)pdfBGwithPops(xdata,BGcurve_interp,output_BG_alldt,frame_time,populations,size_dt,verbose,varargin);
%Set the MLEoptions
mleoptions = statset('mlecustom');
mleoptions.MaxIter = 2000; %Maximum iterations
mleoptions.MaxFunEvals = mleoptions.MaxIter; %Maximum function evaluations
%Perform MLE fitting, and return the parameters
%testing values:
%currently only info: D1
% startparam = [1e-12 15e-9 1 0.00001];
% lowerbound = [1e-13 15e-9 0 0.00];
% upperbound = [1e-10 15e-9 2 1];
startparam = [1e-12 15e-9 1 0.9];
lowerbound = [1e-13 15e-9 0.99 0.0];
upperbound = [1e-10 15e-9 1.01 1];
%[D1 loc_unc 1 dt_1_ratio]


% [output1] = pdfBGwithPops(reconstit_arr,BGcurve_interp,output_BG_alldt,frame_time,populations,size_dt,verbose,{startparam});
% output1(1:5)
% sum(output1)
%

% [parameters] = mle(inputarr, 'pdf',custompdf,'start',startparam,'LowerBound',lowerbound,'UpperBound',upperbound,'Options',mleoptions);
%% Testing two univeriates
clc;
custompdf = @(x,varargin)MLEfitter(x,frame_time,populations,size_dt,verbose,varargin);
inputarr = (reconstit_arr(reconstit_arr<0.5e-6));
start = [0.8e-12];
lb = [1e-15];
ub = [2e-12];
options = statset('MaxIter',500, 'MaxFunEvals',6000);
paramEsts = mle(inputarr, 'pdf',custompdf, 'start',start, ...
    'lower',lb, 'upper',ub, 'options',options)
%%
figure(3);clf(3);
histogram(inputarr);%,'Probability','Normalization');
hold on
xdata = [0:1e-8:2e-6];
D = paramEsts(1);
loc_unc = 15e-9;
frame_time = 0.01;
D=D+loc_unc^2/frame_time;
ydata = xdata.*exp(-(xdata.^2)./(4*D*frame_time));
plot(xdata,ydata*1e8,'k-','LineWidth',4);
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