function [startparam, lowerbound, upperbound] = get_population_parametersFF(populations,paramEsts,HOfitCI,bgratios,start,lb_1pop,ub_1pop,lb_2pop,ub_2pop,performestimationfit)
% Get start param, lower bound, upper bound, based on number of populations
% and estimated params. 
% Required inputs:
% populations: integer - number of populations - 1 to 3
% paramEsts: estimated start params
% HOfitCI: CI output of HistObtain-method, or empty array if not used
% bgratios: start BG ratio values found before
% start: start vals
% lb_1pop: start val for lb_1pop
% ub_1pop: start val for ub_1pop
% lb_2pop: start val for lb_2pop
% ub_2pop: start val for ub_2pop
% performestimationfit: Boolean whether estimation (HO) fit is performed

% Output:
% startparam
% lowerbound
% upperbound

CImultiplier = 50; %additional multiplier on top of 95%CI (i.e. a multiplier of 1 is 95% CI, 2 is ~97.5%, etc.
if populations == 1
    if sum(isnan(HOfitCI)) == 0
        if performestimationfit == 1
            %Based on CI if possible
            startparam = [paramEsts(1), fliplr(bgratios)];
            lowerbound = [min(paramEsts(1)-(paramEsts(1)-HOfitCI(1,1))*CImultiplier,0), fliplr(bgratios)*0.5];
            upperbound = [paramEsts(1)+(HOfitCI(2,1)-paramEsts(1))*CImultiplier fliplr(bgratios)*1.5];
        else %if no estimation fit is used - re-use the ones used for HO, that are the standard ones
            startparam = [start, fliplr(bgratios)];
            lowerbound = [lb_1pop, fliplr(bgratios)*0.8];
            upperbound = [ub_1pop, fliplr(bgratios)*1.2];
        end
    else
        %Cop-out for if there's an issue with CI calculation
        disp('Error in CI calc via first fit - using less accurate lb/ub for second fit')
        startparam = [paramEsts(1), fliplr(bgratios)];
        lowerbound = [paramEsts(1)*0.5, fliplr(bgratios)*0.8];
        upperbound = [paramEsts(1)*1.5, fliplr(bgratios)*1.2];
    end
%     keyboard
    %Set max to BG fractions to just below 1
    upperbound(2:end) = min(upperbound(2:end),1-1e-9);
elseif populations == 2
    if sum(isnan(HOfitCI)) == 0
        if performestimationfit == 1
            %output from earlier: paramEsts and bgratios
            startparam = [paramEsts(1), paramEsts(2), paramEsts(3), fliplr(bgratios)];
            %Based on CI if possible
            lbtripleCI = [max(paramEsts(1)-(paramEsts(1)-HOfitCI(1,1))*CImultiplier,0),...
                max(paramEsts(2)-(paramEsts(2)-HOfitCI(1,2))*CImultiplier,0),...
                max(paramEsts(3)-(paramEsts(3)-HOfitCI(1,3))*CImultiplier,0)];
            ubtripleCI = [paramEsts(1)+(HOfitCI(2,1)-paramEsts(1))*CImultiplier,...
                paramEsts(2)+(HOfitCI(2,2)-paramEsts(2))*CImultiplier,...
                paramEsts(3)+(HOfitCI(2,3)-paramEsts(3))*CImultiplier];
            lowerbound = [lbtripleCI, fliplr(bgratios)*0.5];
            upperbound = [ubtripleCI, fliplr(bgratios)*1.5];
        else %if no estimation fit is used - re-use the ones used for HO, that are the standard ones
            startparam = [start, fliplr(bgratios)];
            lowerbound = [lb_2pop, fliplr(bgratios)*0.5];
            upperbound = [ub_2pop, fliplr(bgratios)*1.5];
        end
    else
        %Cop-out for if there's an issue with CI calculation
        disp('Error in CI calc via first fit - using less accurate lb/ub for second fit')
        lowerbound = [paramEsts(1)*0.5, paramEsts(2)*0.5, paramEsts(3)*0.75, fliplr(bgratios)*0.8];
        upperbound = [paramEsts(1)*1.5, paramEsts(2)*1.5, paramEsts(3)*1.25, fliplr(bgratios)*1.2];
    end
end
end