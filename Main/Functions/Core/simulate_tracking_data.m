%% Function simulate tracking data
%Very basic tracking data simulation
%input parameters:
%Darr, farr: [1-by-4] array of D and fractions (in m2/s)
%num_sim: number of tracks to be simulated
%sim_length: length of each track
%loc_unc: localization uncertainty (in m)
%steptime: time between steps (in s)
%output parameters:
%pos: list x, y, z positions
function pos = simulate_tracking_data(Darr,farr,num_sim,sim_length,loc_unc,steptime)
% JDFittingsimulations

pos = zeros(sim_length*(num_sim),3);
for iter = 1:num_sim %Iterate over number of simulated particles
    %Decide Diff coeff
    randval = rand();
    if randval < farr(1)
        DiffCoeff(iter) = Darr(1);
    elseif randval < sum(farr(1:2))
        DiffCoeff(iter) = Darr(2);
    elseif randval < sum(farr(1:3))
        DiffCoeff(iter) = Darr(3);
    else
        DiffCoeff(iter) = Darr(4);
    end
    for t = 1:sim_length %loop over time points
        %Get new position
        for dim = 1:3
            if t == 1
                pos((iter-1)*(sim_length)+t,dim) = 0;
            else
                pos((iter-1)*(sim_length)+t,dim) = pos((iter-1)*(sim_length)+t-1,dim)+(randn.*sqrt(2*DiffCoeff(iter)*steptime));
            end
        end
    end%end loop over sim length
end %end iteration loop
%Add localizaiton error
localizationerror = normrnd(0,loc_unc, size(pos,1), size(pos,2));
pos = pos + localizationerror; %place where it adds localization uncertainty
end
