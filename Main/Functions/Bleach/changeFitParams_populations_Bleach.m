function [startparam, lowerbound, upperbound, fixRatios2pop] = changeFitParams_populations_Bleach(fitWithBleach, populations, startparam, lowerbound, upperbound, fixRatios2pop, dt_arr)
% Change fitting parameters from [start_BGvalue, bleach_time] to
% [{all_BGvalues}].
% Required inputs:
% fitWithBleach: Boolean whether or not fitted with Bleach
% populations: nr of populations
% startparam
% lowerbound
% upperbound
% fixRatios2pop: Boolean whether the bleach of the 2 populations is fixed
% dt_arr

% Outputs:
% startparam
% lowerbound
% upperbound
% fixRatios2pop - this might be changed

%Change parameters and lower/upper bound in case we fit with bleach
%kinetics rather than fully random BG values
if fitWithBleach
    if populations == 1
        [startparam, lowerbound, upperbound] = changeFitParamsBleach([rand(1)*0.09+0.01,0.0001,5],startparam,lowerbound,upperbound,'1pop');
    elseif populations == 2
        if fixRatios2pop
            [startparam, lowerbound, upperbound] = changeFitParamsBleach([rand(1)*0.03+0.005,0.0001,5],startparam,lowerbound,upperbound,'2pop',[rand(1)*0.05+0.04,0.0001,5]);
            startparam(end) = []; %Fix the 2 population bleach times
            lowerbound(end) = []; %Fix the 2 population bleach times
            upperbound(end) = []; %Fix the 2 population bleach times
        else
            [startparam, lowerbound, upperbound] = changeFitParamsBleach([rand(1)*0.03+0.005,0.0001,5],startparam,lowerbound,upperbound,'2pop',[rand(1)*0.05+0.04,0.0001,5]);
        end
    end
end

% if populations == 2
%     %Parameter initialisation if 2pop ratios are to be fixed or not
%     if fitWithBleach
%         if fixRatios2pop
% 
%         else %not fixing ratios
% 
%         end
%     else %not fitting with Bleach
%         keyboard %figure out
%     end
% %     if fixRatios2pop == 0 && fitWithBleach == 0
% %         %If the ratio for the population is different over DTs, repeat
% %         %the ratio value for the nr of dts
% %         startparam = [startparam(1:2) repmat(startparam(3),1,size(dt_arr,2)) startparam(4:end)];
% %         lowerbound = [lowerbound(1:2) repmat(lowerbound(3),1,size(dt_arr,2)) lowerbound(4:end)];
% %         upperbound = [upperbound(1:2) repmat(upperbound(3),1,size(dt_arr,2)) upperbound(4:end)];
% %     end
% % 
% %     if fitWithBleach
% %         %If fitting with bleach, the ratio over dt will be dependant on the
% %         %bleach time, so we only need to input one ratio
% %         fixRatios2pop = 1;
% %     end
% end
end