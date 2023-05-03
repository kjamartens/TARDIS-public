% % TP/FP: ratio at every dt between red/blue
% % TP+FP: count for every localization at every dt
% % bleach: from curve at dt=1-5
% % density: calc at dt=0 at every frame, average
% Get out TP from TP+FP and TP/FP and density, then see if TP is higher/lower/equal than expected bleach
% Difference of TP <--> bleach is caused by blinking

% parameters = TP/FP, if bleach used, then ca TP/FP from here via
% multiExpFitTARDISBleach
% parameters = bleach, if bleach not used, use TP/FP to get bleach
% (MultiExpFitTARDISBleach for hints)
% JDarrSignalCell = all linkages = TP+FP
% poslist = used for density

function swift_params = obtain_swift_params(parameters,JDarrSignalCell,poslist,maxdist,frame_dist_BG,settingsTARDIS,user_provided_noise_dens,BGarrtotsize,callfromUI)
    dispUIorCommandWindow('Starting obtaining swift parameters....',callfromUI);
    %% Proper density calculation
    user_provided_noise_dens = settingsTARDIS.noiseDensity;
    dens_noise = user_provided_noise_dens; %in loc/um2/frame %
    pos = poslist;
    %Calculte the hypothetical dt1 difference if no blinking is happening
    bleach_time_half = parameters(end)/settingsTARDIS.frame_time; %in frames!
    nr_jumps_per_track = bleach_time_half/log(2)+.5;

    %Get full density out - using 90th percentile on left and right (and
    %top/bottom) to filter out edge-effects
    pc1 = 10;
    pc2 = 100-pc1;
    posxy = pos(pos(:,2)>prctile(pos(:,2),pc1),:);
    posxy = posxy(posxy(:,2)<prctile(posxy(:,2),pc2),:);
    posxy = posxy(posxy(:,3)>prctile(posxy(:,3),pc1),:);
    posxy = posxy(posxy(:,3)<prctile(posxy(:,3),pc2),:);
    %Next var is for the density that is contributed by both the TP and FP
    %(noise) localizations
    dens_signal_and_noise = (size(posxy,1)/(((max(posxy(:,2))-min(posxy(:,2)))*(max(posxy(:,3))-min(posxy(:,3))))*1e12))./max(pos(:,1));
    
    %Get calculated value of spurious noise density out - from full FoV as the
    %noise is spread out over this whole FoV.
    pc1 = 0;
    pc2 = 100-pc1;
    posxy0p = prctile(pos(:,2),pc1);
    posxy0p2 = prctile(pos(:,2),pc2);
    posxy0p3 = prctile(pos(:,3),pc1);
    posxy0p4 = prctile(pos(:,3),pc2);
    %We have the expected density of noise from the user in dens_noise
    %And use that to get a ratio (0-1) of all locs that correspond to noise out
    dens_ratio = (dens_noise*max(pos(:,1))*(((posxy0p2-posxy0p)*(posxy0p4-posxy0p3))*1e12))/size(pos,1);
    
    %We also get the density of only the TP (signal) out.
    dens_signal = dens_signal_and_noise-dens_noise;
    
    %We calculate the number of TP tracks from the TP loc nr and bleach rate
    TPtrack_from_TPloc = length(pos)/nr_jumps_per_track*(1-dens_ratio);
    
    %From this, we calculate the number of linkages at dt1, no blinking assumed, per frame
    dtval = 1;
    TPlinkages_noblink_perframe = TPtrack_from_TPloc*(nr_jumps_per_track-dtval)/max(pos(:,1));
    
    %This is corrected for the TP density, so that we have the number of 
    % linkages at dt1, no blinking assumed, per frame, assuming a density (TP) of 1 loc/um2/frame
    TPlinkages_noblink_perframe_dens1 = TPlinkages_noblink_perframe/dens_signal;
    %%
    %Get the information of how many linkages we have at different dts
    %First of just all the dt-bins that we want
    maxdt = max(settingsTARDIS.dt_arr);
    for i = 1:maxdt
        all_dt_info{i} = JD_Rel(poslist,i,maxdist,0);
        dtbin(i) = i;
    end
    %Then from the frame_dist_BG array (i.e. long frame-delays in TARDIS)
    high_dt_info_length=BGarrtotsize;
    for i = 1:maxdt
        all_dt_info_length(i) = size(all_dt_info{i},1);
    end
    %%
    %Here, we calculate the background localizations at every dt-level,
    %based on the raw data only
    %First, we get the mean localizations per frame
%     zeroFrameJD = JD_Rel(poslist,0,maxdist,0);
    meanlocperframe = size(poslist,1)/max(poslist(:,1));

    %Next, we calculate the ratio of loc combinations inside the search radius vs all locs in
    %frames - needs to be calculated from BG info
    for i = 1:size(frame_dist_BG,2)
        ratio_inSearchRadius_to_allOnFrame_fromBG(i) = (sqrt(high_dt_info_length(i)/(max(poslist(:,1)-frame_dist_BG(i)))))/(meanlocperframe);
    end
    %This is averaged
    ratio_inSearchRadius_to_allOnFrame = mean(ratio_inSearchRadius_to_allOnFrame_fromBG);    
    %Finally, this ratio is used to calculate the expected BG localizations
    %at every dt delay
    for i = 1:maxdt
        calculated_BGlocs(i) = (meanlocperframe*ratio_inSearchRadius_to_allOnFrame)^2*(max(poslist(:,1))-i);
    end
    %We also need the slope and y-intersect for later - currently unused
    %but might be nice?
    BGlocnumber_yintersect = (meanlocperframe*ratio_inSearchRadius_to_allOnFrame)^2*(max(poslist(:,1))-0);
    BGlocnumber_slope = -1*(meanlocperframe*ratio_inSearchRadius_to_allOnFrame)^2;

    % Visualisation here if we ever need it - best to use a dt of >> 3
    % (like 50 or so)
%     keyboard
%     figure(5);clf(5);
%     plot(dtbin,all_dt_info_length/max(poslist(:,1))/dens_signal,'DisplayName','All dt info');
%     hold on; plot([1:maxdt],calculated_BGlocs/max(poslist(:,1))/dens_signal,'k--','DisplayName','Extrapolated BG localizations from high-dt')
%     legend()
%     xlabel('Dt bin')
%     ylabel({'# linkages between localizations per frame,',' corrected for density of TP'})

    dt1_difference_value_gotten = all_dt_info_length(1)/max(poslist(:,1))/dens_signal - calculated_BGlocs(1)/max(poslist(:,1))/dens_signal;

    %%
%     keyboard
    %This value can be compared to the value gotten at the top, and the
    %ratio is the blinking ratio.
    found_blink_ratio = 1-dt1_difference_value_gotten/TPlinkages_noblink_perframe_dens1;
    %%
    % Put in a total array - for now only blinking ratio
    swift_params = max(0,found_blink_ratio);
end 