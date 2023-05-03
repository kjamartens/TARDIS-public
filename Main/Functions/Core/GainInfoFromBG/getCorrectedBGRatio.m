function [SwiftStructure] = getCorrectedBGRatio(JDonlydata,JDarrBG,reconstit_arr,JDarrSignalCell,frame_dist_BG,poslist,settingsTARDIS)

bins = linspace(0,max(JDarrBG),300);

%Get hist counts of signal
histcountsSignal_prob = histcounts(JDonlydata{1}{1},bins,'Normalization','probability');
histcountsSignal_count = histcounts(JDonlydata{1}{1},bins,'Normalization','count');
%Get hist counts of BG
histcountsBG_prob = histcounts(JDarrBG,bins,'Normalization','probability');
histcountsBG_prob = histcountsBG_prob./max(histcountsBG_prob);
histcountsBG_count = histcounts(JDarrBG,bins,'Normalization','count');
%Correct BG for nr of BG frames
histcountsBG_count = histcountsBG_count./size(frame_dist_BG,2);

%Get bins-ratio of this
% histcountsRatio = histcountsBG_count.*histcountsSignal_prob;
histcountsRatio = histcountsBG_prob.*histcountsSignal_count;

Signal_to_BG_ratio_bins = histcountsSignal_prob.*histcountsSignal_count./histcountsBG_count;

%Get loc density
FoVsize = [prctile(poslist(:,2),97.5)-prctile(poslist(:,2),2.5),prctile(poslist(:,3),97.5)-prctile(poslist(:,3),2.5)].*1e6; %in um
% FoVsize = [max(poslist(:,2))-min(poslist(:,2)),max(poslist(:,3))-min(poslist(:,3))].*1e6; %in um
area = FoVsize(1)*FoVsize(2);
loc_density = (size(poslist,1)/max(poslist(:,1)))/area; %in loc/frame/um2

%Calculate corrected value - determined from simulations at 300 bins
baseline = 3.18e-10;

corr_BG_value = max(1,baseline./(mean(Signal_to_BG_ratio_bins).*(1e6*median(JDonlydata{1}{1})^2)*loc_density));

exp_displacement = median(JDonlydata{1}{1})*1e9;
%Get emperically determined blink/bleach:
load('Complexity_fitCurves.mat')
p_blink = ftcurve_blinking(corr_BG_value)/100;
exp_noise_rate = ftcurve_spurious(corr_BG_value)/100;

%Get p_bleach
for dt = 1:size(JDonlydata{1},2)
    size_nr_JDs(dt) = size(JDonlydata{1}{dt},1);
end
size_nr_JDs = size_nr_JDs./size_nr_JDs(1);
bleach_frames = (-1*[1:size(size_nr_JDs,2)]+1)/log2(size_nr_JDs);

%Writing JSON file
%Create a structure with the information
SwiftStructure.exp_displacement = round(exp_displacement,4);
SwiftStructure.exp_noise_rate = round(100*exp_noise_rate,4);
SwiftStructure.p_bleach = round(1/bleach_frames,4);
SwiftStructure.p_blink = round(p_blink,4);
SwiftStructure.p_reappear = round(1-p_blink,4);
SwiftStructure.precision = round(settingsTARDIS.loc_unc*1e9,4);

if settingsTARDIS.StoreSWIFTparameters
    foldernameid = max(strfind(settingsTARDIS.callfromUI.DataLocationEditField.Value,'\'));
    if ~isempty(foldernameid)
        %Get folder and file positions
        folder = settingsTARDIS.callfromUI.DataLocationEditField.Value(1:foldernameid);
        filenameid = max(strfind(settingsTARDIS.callfromUI.DataLocationEditField.Value,'.'));
        filename = [settingsTARDIS.callfromUI.DataLocationEditField.Value(foldernameid+1:filenameid-1) '_SWIFTparams.json'];
        %Write to data
        if ~exist([folder '\TARDIS_Results\'])
            mkdir([folder '\TARDIS_Results\']);
        end
        fileID = fopen([folder '\TARDIS_Results\' filename],'w');
        fprintf(fileID,jsonencode(SwiftStructure,PrettyPrint=true));
        fclose(fileID);
        dispUIorCommandWindow('Stored SWIFT params to JSON file.',settingsTARDIS.callfromUI)
    else
        dispUIorCommandWindow(['Found SWIFT params: p_bleach: ' num2str(round(1/bleach_frames,3)) '; exp_displacement:' num2str(round(exp_displacement,0)) '; exp_noise_rate:'...
        num2str(100*round(exp_noise_rate,4)) '; p_blink:' num2str(round(p_blink,2)) '; p_reappear:' num2str(round(1-p_blink,2)) ';'],settingsTARDIS.callfromUI);
        dispUIorCommandWindow('Error with writing SWIFT values to JSON.',settingsTARDIS.callfromUI)
    end
else
%     dispUIorCommandWindow(['Found SWIFT params: p_bleach: ' num2str(round(1/bleach_frames,3)) '; exp_displacement:' num2str(round(exp_displacement,0)) '; exp_noise_rate:'...
%     num2str(100*round(exp_noise_rate,4)) '; p_blink:' num2str(round(p_blink,2)) '; p_reappear:' num2str(round(1-p_blink,2)) ';'],settingsTARDIS.callfromUI);
end
% res = mean(Signal_to_BG_ratio_bins);
end