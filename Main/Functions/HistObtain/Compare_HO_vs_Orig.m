%% Compare HO (Hist-obtain) vs simulatenous BG-fit method

%% Data loading
% mainFolder = 'C:\Users\Koen Martens\Documents\WUR info - data\Backup H Drive\MATLAB_BIPNAS\RelativeDistance_trackingAlternative\SimulatedTracks\OnePop_anaDDA\';
% dens = {'1over10','1over4','1over2','1','2','3','4','8','10','15','20'};
% Dfree = {'1','2','5','10','22','46','100','215','464','1000'};
% 
% denschoice = 4;
% Dfreechoise = 7; %1,4,7,10
% 
% load([mainFolder 'dens_' dens{denschoice} '_kon_0_koff_0_Dfree_' Dfree{Dfreechoise} '.mat']);

%% Data generation
% dens = 0.05;
nrtrackarr = [100 500 1000 5000];
rep = 20;
densarr = [0.25 1 3];
for nrtrack = 1:size(nrtrackarr,2)
    for ddens = 1:size(densarr,2);
    for rrep = 1:rep
Darr = [1 0 0 0]*1e-12; %Array of Diff Coeff values - always size 4
farr = [1 0 0 0]; %array of fractions of D values - always size 4

num_sim_final = nrtrackarr(nrtrack); %Actual number of simulated particles
num_sim = num_sim_final*4; %number of simulated particles
sim_length = 5; %length of each simulated particle
loc_unc = 15*1e-9;%loc uncertainty in nm
steptime = 10/1000; %steptime in seconds

randomdistributed = 1; %1 = completely random distribution of start points, 0 is very non-random (i.e. structured)
framesize = 5e-6;
tracksstartingperframe = densarr(ddens);%dens*(framesize*1e6)^2/(sim_length); %number of tracks starting per frame

[poslist, ~, ~] = ...
    simulate_tracking_data_with_frame(Darr,farr,num_sim,sim_length,loc_unc,steptime,randomdistributed,framesize,tracksstartingperframe);
            
%% HO method
clear settingsURDA
settingsURDA.frame_dist_BG = [9:1:20];
settingsURDA.startpointBG = 1.5e-6;
settingsURDA.maxdist = 3e-6;
settingsURDA.dt_arr = [1 2 3];
settingsURDA.vis = 1;
settingsURDA.verbose = 0;
settingsURDA.fitvis = 1;
settingsURDA.bgbinningnr = 100;
settingsURDA.linorlogBGsubtract = 'lin';
settingsURDA.linorlogvis = settingsURDA.linorlogBGsubtract;
settingsURDA.minlogpoint = 10^-9;
settingsURDA.frame_time = 0.01;
settingsURDA.loc_unc = 15e-9;
settingsURDA.norm_bins = 10000;
settingsURDA.start = 10.^(-11-rand()*2); %Startpoint of fitting between 1e-11 and 1e-13
settingsURDA.lb = [1e-16]; %LB of fit
settingsURDA.ub = [30e-12]; %UB of fit
settingsURDA.debug = false;

HOfit = URDA_HO_MLE_function(poslist,settingsURDA);

%% Old method
populations = 1;
verbose = 0;
maxdtbins = max(settingsURDA.dt_arr);
startdiff = settingsURDA.start;%(rand()*(20e-12-0.01e-12)+0.01e-12) %start somewhere between 0.01 and 20 um2/s
startparam = [startdiff 15e-6 1 ones(1,maxdtbins)*0.5];
lowerbound = [settingsURDA.lb 1e-6 0.0099 ones(1,maxdtbins)*0];
upperbound = [settingsURDA.ub 50e-6 100.01 ones(1,maxdtbins)*1];
[parameters,t,nll,visualisationinfo] = MLE_BG_subtraction_GeneralPopulations_function(poslist,settingsURDA.frame_dist_BG...
    ,settingsURDA.maxdist,settingsURDA.frame_time,maxdtbins,populations,verbose,startparam,lowerbound,upperbound);

startsettings{nrtrack,ddens,rrep} = settingsURDA;
HOfitcell{nrtrack,ddens,rrep} = HOfit;
origfitcell{nrtrack,ddens,rrep} = parameters;

disp([num2str(nrtrackarr(nrtrack)) ' - ' num2str(densarr(ddens)) ' - ' num2str(rrep) ': HO: ' num2str(HOfit*1e12) '; orig: ' num2str(parameters(1)*1e12)])
    end
    end
end

%% Visualisation of results
clc
figure(19);clf(19);
for i = 1:4
    for j = 1:3
subplot(3,4,(i+(j-1)*4))
for k = 1:20
temparrHO(k) = (HOfitcell{i,j,k});
temparrOrig(k) = (origfitcell{i,j,k}(1));
temparrStartPos(k) = (startsettings{i,j,k}.start);
end
hold on
plot([1:22]-1,ones(1,22)*1e-12,'k-')
plot([1:20],temparrHO,'r*')
plot([1:20],temparrOrig,'b*')
plot([1:20],temparrStartPos,'g*')
axis([0 21 1e-13 1e-11])
set(gca,'YScale','log')
title(['Nr tracks: ' num2str(nrtrackarr(i)) ' Dens: ' num2str(densarr(j))])
    end
end

figure(20);clf(20);

for i = 1:4
    for j = 1:3
subplot(3,4,(i+(j-1)*4))
for k = 1:20
temparrHO(k) = (HOfitcell{i,j,k});
temparrOrig(k) = (origfitcell{i,j,k}(1));
temparrStartPos(k) = (startsettings{i,j,k}.start);
end
%Filter out values that are < 20% off from real value (1e-12)
temparrHO(abs(temparrHO-(1e-12))>((1e-12)*0.2)) = [];
temparrOrig(abs(temparrOrig-(1e-12))>((1e-12)*0.2)) = [];
hold on
ls = linspace(0.9e-12,1.1e-12,10);
h=histogram(temparrHO,ls);
h.FaceColor = 'r';
h.EdgeColor = 'none';
h.FaceAlpha = 0.4;
h=histogram(temparrOrig,ls);
h.FaceColor = 'b';
h.EdgeColor = 'none';
h.FaceAlpha = 0.4;
% disp(['Nr tracks: ' num2str(nrtrackarr(i)) ' Dens: ' num2str(densarr(j)) ' STD HO: ' num2str(std(temparrHO)) ' STD Orig: ' num2str(std(temparrOrig))])
disp(['Nr tracks: ' num2str(nrtrackarr(i)) 9 ' Dens: ' num2str(densarr(j)) 9 'MeanErr HO: ' num2str(mean(abs(temparrHO-(1e-12)))*1e12*100) 9 ' MeanErr Orig: ' num2str(mean(abs(temparrOrig-(1e-12)))*1e12*100) 9 ' Ratio orig/HO: ' num2str(mean(abs(temparrOrig-(1e-12)))/mean(abs(temparrHO-(1e-12))))])
% plot([1:22]-1,ones(1,22)*1e-12,'k-')
% plot([1:size(temparrHO,2)],temparrHO,'r*')
% plot([1:size(temparrOrig,2)],temparrOrig,'b*')
% plot([1:20],temparrStartPos,'g*')
% axis([0 21 1e-13 1e-11])
% set(gca,'YScale','log')
title(['Nr tracks: ' num2str(nrtrackarr(i)) ' Dens: ' num2str(densarr(j))])
    end
end