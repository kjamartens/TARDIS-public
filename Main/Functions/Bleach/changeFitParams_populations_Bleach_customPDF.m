function [startparam, lowerbound, upperbound] = changeFitParams_populations_Bleach_customPDF(fitWithBleach, nrparamskeep, startparam, lowerbound, upperbound)
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
    s = [0 startparam(nrparamskeep+1:end)];
    l = [0 lowerbound(nrparamskeep+1:end)];
    u = [0 upperbound(nrparamskeep+1:end)];
    [ss, ll, uu] = changeFitParamsBleach([rand(1)*0.09+0.01,0.0001,5],s,l,u,'1pop');

    startparam = [startparam(1:nrparamskeep) ss(2:end)];
    lowerbound = [lowerbound(1:nrparamskeep) ll(2:end)];
    upperbound = [upperbound(1:nrparamskeep) uu(2:end)];
end

end