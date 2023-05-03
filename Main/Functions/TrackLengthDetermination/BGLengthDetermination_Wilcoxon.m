function BGframe = BGLengthDetermination_Wilcoxon(pos,nrsamplepoints,nriterationsformedian,maxdist,visualisation)
% Input:
%   pos: position data - columns should be frame-x-y. Units in frame or m.
%
%   nrsamplepoints: the number of points the algorithm should sample at.
%   More points = better statistics, but slower. Recommended at 16 or 25.
%
%   nriterationsformedian: the number of iterations after which a median
%   is performed for determining the frame at which only BG is there.
%   Recommended 10-ish.
%   
%   maxdist: maximum distance to probe at. Probably use the same value as
%   used in TARDIS. Units in m.
%
%   visualisation: set to 1 for figures
% Output:
%   BGframe = lowest determined frame at which it's only Background.

tic
%We hard-code some min and max deltaTs.
minDeltaT = 10; %Proper assumption that there will always be tracks longer than 10 frames,
%and using too low minDeltaT removes effective info (i.e. logspace from
%1-20ish is worthless but uses 6 bins or so).
maxDeltaT = max(pos(:,1)/4); %Quarter of total nr of frames is more or less the 
%safest position at which we're sure all data is real BG data

%We sample between the min and max deltaT logarithmically
deltatarr = round(logspace(log10(minDeltaT),log10(maxDeltaT),nrsamplepoints));
%%
%Create indexing for frame positions
aarr = cumsum(accumarray(pos(:,1),1));
%Loop over all deltaTs that we want to look at
for i = 1:size(deltatarr,2)
    %Start a counter
    bigcounter = 1;
    %Start an empty array
    fulldist = zeros(10000000,1);
    reps = 1; %This value is determining to look at only e.g. frame 1-51, or also 1-52, 1-53 (at reps = 3). I think reps = 1 is just good
    for rrep = 1:reps
        %Loop over all frames
        %Skip frame 1, and end early enough
        for f = 2:length(aarr)-1-deltatarr(i)-reps
            %Get the information of these two frames - this basically pulls
            %out data of frame1 vs frame2, assuming that's f and
            %f+deltaarr(i).
            frameinfo1 = pos(aarr(f-1)+1:aarr(f),:);
            frameinfo2 = pos(aarr(f-1+deltatarr(i))+1:aarr(f+deltatarr(i)),:)+rrep-1;

            %Get distance between all locs in these frames
            %For this, we create an empty array or correct size
            tempdist = zeros(size(frameinfo1,1)*size(frameinfo2,1),1);
            %And create a counter
            counter = 1;
            %Loop over all locs in the first frame
            for l1 = 1:size(frameinfo1,1)
                %Loop over all locs in the second frame
                for l2 = 1:size(frameinfo2,1)
                    %Calculate the distance between the locs and store it
                    %if it's small enough
                    dist = sqrt((frameinfo1(l1,2)-frameinfo2(l2,2))^2+(frameinfo1(l1,3)-frameinfo2(l2,3))^2);
                    if dist < maxdist
                        tempdist(counter) = dist;
                        counter = counter+1;
                    end
                end
            end
            %Fill the large arrays with this distance array
            fulldist(bigcounter:bigcounter+size(tempdist,1)-1) = tempdist;
            %Increase the required counter
            bigcounter = bigcounter+size(tempdist,1)+1;
        end
    end
    %Remove zero-entries (original array is created too large)
    fulldist(fulldist==0) = [];
    %Store as Cell
    fulldistcfull{i} = fulldist;
end
%% Random sampling and Wilcoxon test
%Next, we take random samples from these datasets, based on the smallest
%size.
%First we determine the smallest length of these arrays, and store it in
%randsamp
randsamp = 9e10;
for k = 1:nrsamplepoints
    if size(fulldistcfull{k},1) < randsamp
        randsamp = size(fulldistcfull{k},1);
    end
end

%Now we start determining if these populations are different than the BG
%population.
%We iterate over the number of iterations inputted.
for it = 1:nriterationsformedian
    %We create the randomly-sampled datasets
    for c = 1:size(deltatarr,2)
        fulldistc{c} = datasample(fulldistcfull{c},randsamp,'Replace',false);
    end
    %We keep track on first time the Wilcoxon test fails
    firstinstanceWilcoxontestfail(it) = 0;
    %We perform the W-test at all datasets (even last one, which will never
    %fail, since it's the same dataset...)
    for c = 1:size(deltatarr,2)
        %Wilcoxon test is in 'ranksum'. We test if the lower-deltaT dataset
        %has a left-shifted (i.e. lower values) median value than the
        %dataset at very high deltaT. I believe it should always be
        %left-shifted, but might be changed at some point.
        %It only passes with a p = 0.01 value.
        [p(c),h] = ranksum(fulldistc{c},fulldistc{size(deltatarr,2)},'Alpha',1e-2,'tail','left');
        %We keep track if this is the first one that fails (i.e. showcasing
        %that it's the 'same' distribution statistically)
        if h == 0
            if firstinstanceWilcoxontestfail(it) == 0
                firstinstanceWilcoxontestfail(it) = c;
            end
        end
    end
    %We can display this if wanted - turned off for now
%     disp(deltatarr(firstinstanceWilcoxontestfail(it)))
end

%Now we take the median value of when it first fails, and that's our BG
%start point
BGframe = deltatarr(floor(median(firstinstanceWilcoxontestfail)));
%We can display if needed - turned off for now
% disp(['-- ' num2str(BGframe) ' --'])
toc;
%%
%Visualisation if we want. Note that it's slow, because I perform the
%Wilcoxon test again everywhere...
if visualisation
    figure(3);clf(3);
    bininfo = linspace(0,maxdist,20);
    for i = 1:size(deltatarr,2)
        subplot(ceil(sqrt(nrsamplepoints)),ceil(sqrt(nrsamplepoints)),i)
        histogram(fulldistc{i},bininfo,'Normalization','probability','FaceColor',[0.6 0.6 0.6])
        hold on
        histogram(fulldistc{size(deltatarr,2)},bininfo,'DisplayStyle','stairs','Normalization','probability','LineWidth',2,'EdgeColor','k')
        
        [p,h] = ranksum(fulldistc{i},fulldistc{size(deltatarr,2)},'Alpha',1e-2,'tail','left');
        failorsucceedtxt = 'Fail';
        if h == 1
            failorsucceedtxt = 'Succeed';
        end
        title(['\Deltat: ' num2str(deltatarr(i)) ' frames: ' failorsucceedtxt])
    %     title(['\Deltat: ' num2str(deltatarr(i)) ' frames'])
        ylabel('pdf')
        xlabel('JD (nm)')
    end
    %
    
    figure(4);clf(4);
    bininfo = linspace(0,2e3,50);
    for i = 1:size(deltatarr,2)
        subplot(ceil(sqrt(nrsamplepoints)),ceil(sqrt(nrsamplepoints)),i)
        ll = cdfplot(fulldistc{i});
        ll.Color = [1 0 0];
        ll.LineStyle = '-';
        ll.LineWidth = 2;
        hold on
        ll = cdfplot(fulldistc{size(deltatarr,2)});
        ll.Color = [0 0 0];
        ll.LineStyle = '--';
        ll.LineWidth = 2;
    %     [h,p(i)] = kstest2(fulldistc{i},fulldistc{size(deltatarr,2)},'Alpha',1e-2,'Tail','larger');
        [p,h] = ranksum(fulldistc{i},fulldistc{size(deltatarr,2)},'Alpha',1e-2,'tail','left');
        failorsucceedtxt = 'Fail';
        if h == 1
            failorsucceedtxt = 'Succeed';
        end
        title(['\Deltat: ' num2str(deltatarr(i)) ' frames: ' failorsucceedtxt])
    end
end
end