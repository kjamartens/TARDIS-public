clear all
%Load data
% filename = '\\IFMB-NAS\AG Endesfelder\Data\Koen\Articles\2017_miCube\DataNatComm\pNonTarget_day1_locData_ThSTORMinput.csv';
% locdata = commonFunctions.readCSVThunderSTORM(filename);

% load('\\IFMB-NAS\AG Endesfelder\Data\Koen\Articles\2020_relativeDisplacementTracking\ManuscriptFigureData\anaDDANoise\aDDANoise_20_20_2_oLoc_0_sLoc_0.mat')
% pos = commonFunctions.readCSVThunderSTORM('\\IFMB-NAS\AG Endesfelder\Data\Koen\transfer\2mW_nonFTMThSTORM.csv');
% locdata = pos;
% locdata(:,2:3) = locdata(:,2:3)*1e-9;


load('\\IFMB-NAS\AG Endesfelder\Data\Koen\Scientific\MATLAB\GITHUB\LongLivedPhotophysics\SimulatedTracks\BleachTracks_noBlinkingIncorp_02sechalfBleachTime_5ktracks.mat')

% load('\\IFMB-NAS\AG Endesfelder\Data\Koen\Scientific\MATLAB\GITHUB\LongLivedPhotophysics\SimulatedTracks\BleachTracks_noBlinkingIncorp_05sechalfBleachTime_5ktracks.mat')
locdata = pos;
tic
% deltatarr = [1 5 10 50 100 500 1000 5000 10000];
dv = 3;
maxv = 500;
% deltatarr = round(logspace(log10(1),log10(max(locdata(:,1)/4)),dv^2));
deltatarr = round(logspace(log10(1),log10(maxv),dv^2));
reps = 1;
maxdist = 2e-6;
%%
for i = 1:size(deltatarr,2)
    %Loop over all frames
    aarr = cumsum(accumarray(locdata(:,1),1));
    bigcounter = 1;
    fulldist = zeros(10000000,1);
    for rrep = 1:reps
%         for f = 2:length(aarr)-1-max(deltatarr)-reps
        for f = 2:length(aarr)-1-deltatarr(i)-reps
            frameinfo1 = locdata(aarr(f-1)+1:aarr(f),:);
            frameinfo2 = locdata(aarr(f-1+deltatarr(i))+1:aarr(f+deltatarr(i)),:)+rrep-1;
            %Get distance between all locs
            tempdist = zeros(size(frameinfo1,1)*size(frameinfo2,1),1);
            counter = 1;
            for l1 = 1:size(frameinfo1,1)
                for l2 = 1:size(frameinfo2,1)
                    dist = sqrt((frameinfo1(l1,2)-frameinfo2(l2,2))^2+(frameinfo1(l1,3)-frameinfo2(l2,3))^2);
                    if dist < maxdist
                    tempdist(counter) = dist;
                    counter = counter+1;
                    end
                end
            end
            fulldist(bigcounter:bigcounter+size(tempdist,1)-1) = tempdist;
            bigcounter = bigcounter+size(tempdist,1)+1;
        end
    end
    fulldist(fulldist==0) = [];
    fulldistcfull{i} = fulldist;
end
%%
randsamp = round(size(fulldistcfull{dv^2},1)/1);%
for it = 1:10
    for c = 1:size(deltatarr,2)
        fulldistc{c} = datasample(fulldistcfull{c},randsamp,'Replace',false);
    end
    firstinstanceKStestfail(it) = 0;
    for c = 1:size(deltatarr,2)
%         [h,p(c)] = kstest2(fulldistc{c},fulldistc{size(deltatarr,2)},'Alpha',1e-2,'Tail','larger');
        [p(c),h] = ranksum(fulldistc{c},fulldistc{size(deltatarr,2)},'Alpha',1e-2,'tail','left');
        if h == 0
            if firstinstanceKStestfail(it) == 0
                firstinstanceKStestfail(it) = c;
            end
        end
    end
    disp(deltatarr(firstinstanceKStestfail(it)))
end

disp(['-- ' num2str(deltatarr(floor(median(firstinstanceKStestfail)))) ' --'])
toc
%%
figure(3);clf(3);
bininfo = linspace(0,maxdist,20);
for i = 1:size(deltatarr,2)
    subplot(dv,dv,i)
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
    subplot(dv,dv,i)
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
%     title([failorsucceedtxt ' - ' num2str(hW) ' - ' num2str(pW)])
end


% figure(5);clf(5);
% plot(deltatarr,p)
% set(gca,'XScale','log')
% set(gca,'YScale','log')