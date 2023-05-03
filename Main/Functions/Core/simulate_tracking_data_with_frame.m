%% Function simulate_tracking_data_with_frame
%Very basic tracking data simulation
%input parameters:
%Darr, farr: [1-by-4] array of D and fractions (in m2/s)
%num_sim: number of tracks to be simulated
%sim_length: length of each track
%loc_unc: localization uncertainty (in m)
%steptime: time between steps (in s)
%randomdistributed: Boolean value: TRUE for fully random, FALSE for
%non-random (specified in this script)
%ImageSize: Size of image (in m) - midpoint ~ 5e-6
%tracksstartingperframe: number of tracks starting every frame
%output parameters:
%poswithframeordered: list of frame, x, y, z positions, ordered by frame
%poswithframe: list of frame, x, y, z positions, disordered by frame
%pos: list x, y, z positions

%Scripts dependancy: simulate_tracking_data

function [poswithframeordered, poswithframe, pos] = simulate_tracking_data_with_frame...
    (Darr,farr,num_sim,sim_length,loc_unc,steptime,randomdistributed,ImageSize,tracksstartingperframe)

pos=simulate_tracking_data(Darr,farr,num_sim,sim_length,loc_unc,steptime);
trackstartposIDs =  [1:sim_length:(num_sim*sim_length)-1];
%Now add random starting point to every track
for i = trackstartposIDs
    if randomdistributed
        randval = rand(1,3);
        pos(i:i+sim_length-1,:) = pos(i:i+sim_length-1,:) + randval.*ImageSize;
    else
        %Stupid test to see if non-random distribution does something now
        done = 0;
        while done == 0
            randval = rand(1,3);
%             if (randval(1)< 0.2) || (randval(1) > 0.8)
            if (mod(randval(1),0.2)<(0.1*(randval(1)))) && (mod(randval(2),0.2)<0.05)
                pos(i:i+sim_length-1,:) = pos(i:i+sim_length-1,:) + randval.*ImageSize;
                done = 1;
            end
        end
    end
end

%Now add random frame start to every track
maxFrames = num_sim/tracksstartingperframe; %10 tracks starting every frame
poswithframe = [zeros(size(pos,1),1) pos];
% keyboard
for i = trackstartposIDs
    randval = ceil(rand*maxFrames);% ceil(i*1.5);%
    for j = 1:sim_length
        poswithframe(i+j-1,1) = randval+j-1;
    end
end

%Sort by frame
[~,sortorder] = sort(poswithframe(:,1));
poswithframeordered = poswithframe(sortorder,:);
end