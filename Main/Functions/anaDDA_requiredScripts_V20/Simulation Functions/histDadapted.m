function Dfit = histDadapted(tracks, dT, Numpop, noplot)
% function D = histD(tracks, params) calculates the diffusion coefficient
% for each track in "tracks" and plots a histogram of D.
%
% Stephan Uphoff. 09.09.11
%
% save D histogram in text file- to combine several data later on
% AP 30.01.12
%
% specify textfile of Dcoeff: number of localized molecules, number of tracked molecules, parameter settings
% AP 03.06.12
keyboard
pixel = 1; % length per pixel
%dT = 0.01; % time per frame
sigma = 0.04; % localization noise
sig = sigma^2/dT;
DhistMinSteps = 4; % minimum number of steps for a track to be analyzed
plotfit = 1; % determines whether you also want to include fit of histogram to theoretical distributions 
DhistMaxSteps = 200;
rangeD = 0:0.05:5; % D range for the histogram
nMolecules = max(tracks(:,3));
%noplot = true;

%% if there is cell data present, initialise some variables for generating extra histograms.
cellData=false;

%Option
MSD = zeros(nMolecules,1);
MSDfit = zeros(10*nMolecules,1);
kk = 1;
mm = 1;


%% Part added by JV to enable live video addition of Histogram to video 
table = tabulate(tracks(:,4));
selectedmolecules = table(table(:,2)>DhistMinSteps & table(:,2)<DhistMaxSteps);
MSDold = 0;
threshold = 0.1;
for ii = 1:numel(selectedmolecules)
    xx = find(tracks(:,4)==selectedmolecules(ii));
    l = 0;
    for jj = 2:numel(xx)-1   
        displacement = ((tracks(xx(jj+1),1) - mean(tracks(xx(1:jj),1)))^2 +...
                    (tracks(xx(jj+1),2) - mean(tracks(xx(1:jj),2)))^2);
            if sqrt(displacement) > threshold
            try 
           displacementnextframe = ((tracks(xx(jj+2),1) - mean(tracks(xx(1:jj),1)))^2 +...
                    (tracks(xx(jj+2),2) - mean(tracks(xx(1:jj),2)))^2);    
                if sqrt(displacementnextframe) > threshold
                    break
                else
                 l = l+1;
                end
            catch
            continue
             end
            else
            l = l+1;
            end
    end
     boundtime(ii) = l;
end
for i = 1:15
number = find(boundtime>(i-1));     
numberboundtime(i) = numel(number);
end

%% Start of algorithm 
for ii = 1:nMolecules
    xx = find(tracks(:,4)==ii);
    Amount(ii) = numel(xx);
    

%%   Option1: Average D over all MSD of track (used all the time)
    if numel(xx) > DhistMinSteps && numel(xx) <DhistMaxSteps %
% sum all squared displacement in the track
    for jj = 1:numel(xx)-1
            if tracks(xx(jj+1),3)-tracks(xx(jj),3)==2
                    MSD(kk) = MSD(kk) + (((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                        (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2)/2);
                    else   
                        MSD(kk) = MSD(kk) + ((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                        (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);
                    end   
            end   
    MSD(kk) = MSD(kk)/jj; % mean square displacement   

    if plotfit == 1
% Algorithm used to check if there is a consecutive series somewhere in
% this series
t =  tracks(xx(1:Amount(ii)),3)'; 
N = DhistMinSteps; % Required number of consecutive numbers following a first one
 x = diff(t)==1;
 f = find([false,x]~=[x,false]);
 g = find(f(2:2:end)-f(1:2:end-1)>=N,1,'first');
 first_t = t(f(2*g-1)); % First t followed by >=N consecutive numbers
 if isempty(first_t)
     if sum(x) >= DhistMinSteps
        range2 = find(x==1);
        range2 = range2(1:DhistMinSteps);
     else
        range2 = 1:DhistMinSteps;
     end
 else 
start = find(tracks(xx,3) == first_t);
range2 = start:(start+DhistMinSteps-1);
end
% sum all squared displacement in the track
 
    for jj = range2     
                    if tracks(xx(jj+1),3)-tracks(xx(jj),3)==2
                    MSDfit(mm) = MSDfit(mm) + (((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                        (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2)/2);
                    else   
                        MSDfit(mm) = MSDfit(mm) + ((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
                        (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);
                    end
       end
       mm = mm+1;
    end
    end     
%% if celldata is present, log the cell in which it is located.
       if cellData
            ff=find(tracks(:,4)==ii,1,'First');
            cell=tracks(ff,5);
            TracksPCell(cell) = TracksPCell(cell)+1;
            CellTrack(kk)=cell;
        end
        
    kk = kk + 1;
    
end
a = find(MSD>0);
MSD = MSD(a);
% %   Option2: D for all MSDs (if D(time) needed only)
%     if numel(xx) > DhistMinSteps
%         MSD = zeros(numel(xx),3);
%         % calculate step wise all MSDs and apparent D coefficients
%         for jj = 1:numel(xx)-1
%             
%             MSD(jj,1) = ((tracks(xx(jj+1),1) - tracks(xx(jj),1))^2 +...
%                 (tracks(xx(jj+1),2) - tracks(xx(jj),2))^2);
%             MSD(jj,1) = MSD(jj,1) * pixel^2;
%             MSD(jj,2) = MSD(jj,1)/(4*dT)-sigmaNoise^2*pixel^2/dT; %diffusion coeff with localization noise correction
%            
%             MSD(jj,3) = MSD(jj,1)/(4*dT);
%         end
%         MeanSD = [MeanSD MSD(jj,2)];
%         trackLength = numel(xx);
%         trackID = tracks(xx(1),4);
%     end
 MSDfit = MSDfit/DhistMinSteps;
%% JV added to calculate distribution of track sizes
 range = 0:200;
 NumMolecules = numel(find(MSD>0));
 aveAmount = mean(Amount);
 Amount(Amount==1) = [];
 aveAmountmorethanone = mean(Amount);
 maxAmount = max(Amount);
 %figure;
 Amount = Amount -1;
 %Tracklength = hist(Amount, range);
% xlabel('Track length')
% ylabel('Count')
% axis([-1 20 0 inf]);
%xlswrite('tracklength',Tracklength) %Writes Excel file for the tracklength

%%
if cellData
figure;
hist(TracksPCell, range)
xlabel('Tracks per cell')
ylabel('Count')
axis([-1 51 0 inf]);
end

MSD(kk:end) = []; % delete unused rows

if cellData
CellTrack(kk:end) = []; % delete unused rows
Cells = unique(CellTrack); % find all cell numbers that contain tracks

end

MSD = MSD * pixel^2; % convert from pixel to length units
MSDfit = MSDfit*pixel^2;
% MSDfit = MSD;
% calculate D from MSD and correct for localization noise
%Option1
D = MSD/(4*dT); %- sigmaNoise^2*pixel^2/dT;
Dfit = MSDfit/(4*dT); %- sigmaNoise^2*pixel^2/dT;
Dfit = Dfit(find(Dfit>0));
%Option2
%D = mean(MSD(:,2)); %only for plotting purpose

%Option2
%savename = ['AppDiffCoeff-TrackLength' num2str(trackLength) '-TrackID' num2str(trackID) '_2012-11-20' '.mat'];
%save(savename,'MSD');

if cellData
    % determine the average diffusion coefficient per cell
    for i=1:size(Cells)
        index = CellTrack==Cells(i);
        D_S = D(index);
        AVG_D_S = mean(D_S);
        D_Cell(Cells(i))=AVG_D_S;
    end
end
% plot histogram of single-molecule diffusion coefficients over rangeD
if noplot == false
figure;
logD=real(log10(D));
rangelogD = -3:0.1:1;
hist(logD, rangelogD)
xlabel('Log diffusion coefficient [um^2/s]');
ylabel('Histogram count');
figure;
histogram(D,rangeD)
xlabel('diffusion coefficient [um^2/s]');
ylabel('Histogram count');

if cellData
    % plot the additional histograms that contain cell-based data.
figure;
range2=1:nCells;
plot(range2,D_Cell)
xlabel('Cell number');
ylabel('Average diffusion coefficient [um^2/s]');

figure;
plot(TracksPCell,D_Cell,'.')
xlabel('Track number in the cell');
ylabel('Average diffusion coefficient in the cell [um^2/s]');

figure;
plot(cellAreas,D_Cell,'.')
xlabel('Cell size (area) [um^2]');
ylabel('Average diffusion coefficient [um^2/s]');

figure;
plot(cellVolumes,D_Cell,'.')
xlabel('Cell size (volume) [um^3]');
ylabel('Average diffusion coefficient [um^2/s]');
CombinedCellData = num2cell([range2' TracksPCell D_Cell cellAreas cellVolumes]);
Titles = {'Cell number' 'Number of Tracks' 'Average Diffusion Coefficient' 'Cell Area' 'Cell Volume'};
CombinedCellData = [Titles;CombinedCellData];
savename = [params.tracksPathname params.tracksFilename '_CellData_steplength' num2str(DhistMinSteps(1)) '.xlsx'];
keyboard
xlswrite(savename, CombinedCellData);
end

% modified by SVDW:
% made this optional because it takes a long time.
%if writeXLS
Histvalues = hist(D,rangeD);
Histvalueslog = hist(logD, rangelogD);

%% JV added 27-10-2017 to allow for fitting of the histogram with probability Distributions (Stracy et al. 2015)
%This code also depends on partitioning the tracks into certain lengths.
%That's why it works with a separate variable called Dfit that was
%partitioned earlier in the code. If you don't want to run the script make
%plotfit into different value than 1. 
if plotfit == 1
    data = Dfit;
    x = [0:0.01:5];
    rangeDdifference = rangeD(2)-rangeD(1);
    custompdfs{1} = ') (((DhistMinSteps/lambda).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda))/factorial(DhistMinSteps-1))';
    custompdfs{2} = ') c1.*(((DhistMinSteps/lambda1).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda1))/factorial(DhistMinSteps-1))+(1-c1).*(((DhistMinSteps/lambda2).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda2))/factorial(DhistMinSteps-1))';
    custompdfs{3} = ') c1.*(((DhistMinSteps/lambda1).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda1))/factorial(DhistMinSteps-1))+ c2.*(((DhistMinSteps/lambda2).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda2))/factorial(DhistMinSteps-1))+(1-c1-c2).*(((DhistMinSteps/lambda3).^DhistMinSteps).*(data.^((DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda3))/factorial(DhistMinSteps-1)))';
    custompdfs{4} = ') c1.*(((DhistMinSteps/lambda1).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda1))/factorial(DhistMinSteps-1))+ c2.*(((DhistMinSteps/lambda2).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda2))/factorial(DhistMinSteps-1))+c3.*(((DhistMinSteps/lambda3).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda3))/factorial(DhistMinSteps-1))+(1-c1-c2-c3).*(((DhistMinSteps/lambda4).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda4))/factorial(DhistMinSteps-1))';
    ct{1} = ',c1';
    ct{2} = ',c2';
    ct{3} = ',c3';
    lam{1}= ',lambda1';
    lam{2} = ',lambda2';
    lam{3} = ',lambda3';
    lam{4} = ',lambda4';
    
         lambda1 = NaN;
         lambda2 = NaN;
         lambda3 = NaN;
         lambda4 = NaN;
         c1 = NaN;
         c2 = NaN;
         c3 = NaN;
    
    lambda = [lambda1 lambda2 lambda3 lambda4];
    c = [c1 c2 c3];
    lamx = find(isnan(lambda));
    lamxinv = find(~isnan(lambda));
    lamx = lamx(lamx<=Numpop);
    lamxinv = lamxinv(lamxinv<=Numpop);
    cx = find(isnan(c));
    cxinv = find(~isnan(c));
    cx = cx(cx<Numpop);
    cxinv = cxinv(cxinv<Numpop);
   
    funcstring = strcat('@(data', ct{cx},lam{lamx},custompdfs{Numpop});
    custompdf = eval(funcstring);
    cstart = [0.4,0.5,0.3];
    lamstart = [0.1, 0.5, 1, 2.6];
    totstart = horzcat(cstart(cx), lamstart(lamx));
    if Numpop == 1
        custompdf = @(data,lambda) (((DhistMinSteps/lambda).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda))/factorial(DhistMinSteps-1));
    end
    try
    [Param, CI]= mle(data,'pdf', custompdf, 'start', totstart);
    catch
    keyboard
     totstart = horzcat(cstart(cx), lamstart(lamx));
     [Param, CI]= mle(data,'pdf', custompdf, 'start', totstart);
    end
    
    Paramtemp = zeros(2*Numpop,1);
    Paramtemp(cx) = Param(1:numel(cx));
    Paramtemp(cxinv) = c(cxinv);
    Paramtemp(Numpop) = 1-sum(Paramtemp(1:Numpop-1));
    Paramtemp((lamx+Numpop)) = Param(numel(cx)+1:numel(cx)+numel(lamx));
    Paramtemp(lamxinv+Numpop) = lambda(lamxinv);
    Param = Paramtemp
    figure();
    hold on
    singlecustompdf = @(data,lambda) (((DhistMinSteps/lambda).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda))/factorial(DhistMinSteps-1));
    if Numpop == 1
    custompdf = @(data,lambda) (((DhistMinSteps/lambda).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda))/factorial(DhistMinSteps-1));
    histogram(Dfit,rangeD, 'Normalization','probability', 'Facecolor', [0.6 0.6 0.6], 'Edgecolor',[0.6 0.6 0.6] )
    plot(x,rangeDdifference*custompdf(x,Param(2)),'LineWidth',2, 'Color','black' )
    elseif Numpop == 2
    custompdf = @(data,c1,lambda1,lambda2) c1.*(((DhistMinSteps/lambda1).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda1))/factorial(DhistMinSteps-1))+(1-c1).*(((DhistMinSteps/lambda2).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda2))/factorial(DhistMinSteps-1));   
    histogram(Dfit,rangeD, 'Normalization','probability', 'Facecolor', [0.6 0.6 0.6], 'Edgecolor',[0.6 0.6 0.6] )
    plot(x,rangeDdifference*custompdf(x,Param(1), Param(3), Param(4)),'LineWidth',2, 'Color','black' )
    plot(x,rangeDdifference*Param(1)*singlecustompdf(x,Param(3)),'--', 'LineWidth',2.5, 'Color', [0.3 0.3 0.8]);
    plot(x,rangeDdifference*Param(2)*singlecustompdf(x,Param(4)),'--','LineWidth',2.5, 'Color', [0.3 0.8 0.3]);
    elseif Numpop == 3
    custompdf = @(data, c1, c2, lambda1, lambda2, lambda3) c1.*(((DhistMinSteps/lambda1).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda1))/factorial(DhistMinSteps-1))+ c2.*(((DhistMinSteps/lambda2).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda2))/factorial(DhistMinSteps-1))+(1-c1-c2).*(((DhistMinSteps/lambda3).^DhistMinSteps).*(data.^((DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda3))/factorial(DhistMinSteps-1)));    
    histogram(Dfit,rangeD, 'Normalization','probability', 'Facecolor', [0.6 0.6 0.6], 'Edgecolor',[0.6 0.6 0.6] )
    plot(x,rangeDdifference*custompdf(x,Param(1), Param(2), Param(4),Param(5), Param(6)),'LineWidth',2, 'Color','black' )
    plot(x,rangeDdifference*Param(1)*singlecustompdf(x,Param(4)),'--', 'LineWidth',2.5, 'Color', [0.3 0.3 0.8]);
    plot(x,rangeDdifference*Param(2)*singlecustompdf(x,Param(5)),'--','LineWidth',2.5, 'Color', [0.3 0.8 0.3]);
    plot(x,rangeDdifference*Param(3)*singlecustompdf(x,Param(6)),'--','LineWidth',2.5, 'Color', [0.8 0.3 0.3]);
        
    elseif Numpop == 4    
    custompdf = @(data,c1, c2, c3, lambda1, lambda2, lambda3, lambda4) c1.*(((DhistMinSteps/lambda1).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda1))/factorial(DhistMinSteps-1))+ c2.*(((DhistMinSteps/lambda2).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda2))/factorial(DhistMinSteps-1))+c3.*(((DhistMinSteps/lambda3).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda3))/factorial(DhistMinSteps-1))+(1-c1-c2-c3).*(((DhistMinSteps/lambda4).^DhistMinSteps).*(data.^(DhistMinSteps-1)).*exp(-(DhistMinSteps.*data/lambda4))/factorial(DhistMinSteps-1));
    histogram(Dfit,rangeD, 'Normalization','probability', 'Facecolor', [0.6 0.6 0.6], 'Edgecolor',[0.6 0.6 0.6] )
    plot(x,rangeDdifference*custompdf(x,Param(1), Param(2), Param(3), Param(5), Param(6), Param(7), Param(8)),'LineWidth',2,  'Color','black' )
    plot(x,rangeDdifference*Param(1)*singlecustompdf(x,Param(5)),'--', 'LineWidth',2.5, 'Color', [0.3 0.3 0.8]);
    plot(x,rangeDdifference*Param(2)*singlecustompdf(x,Param(6)),'--','LineWidth',2.5, 'Color', [0.3 0.8 0.3]);
    plot(x,rangeDdifference*Param(3)*singlecustompdf(x,Param(7)),'--','LineWidth',2.5, 'Color', [0.8 0.3 0.3]);
    plot(x,rangeDdifference*Param(4)*singlecustompdf(x,Param(8)),'--','LineWidth',2.5, 'Color', [0.5 0.5 0]);
    end    
    end
        xlabel('Diffusion Coefficient \mum^2/s')
        ylabel('Fraction of tracks')
        legend({'Data', 'Fit'},'FontSize',14)
        set(gca,'fontsize',14)
        hold off   
AmountofMolecules(3,1) = numel(Dfit)
averageD = mean(Dfit)
%filename = [params.tracksPathname params.tracksFilename '_Dhist' '.xlsx'];
%filenamelog = [params.tracksPathname params.tracksFilename '_Dhistlog' '.xlsx'];
%filenamerestdata = [params.tracksPathname params.tracksFilename '_Parameters' '.xlsx']; 
%Matrx = horzcat(rangeD', Histvalues');
%Matrxlog = horzcat(rangelogD', Histvalueslog');
%xlswrite(filename, Matrx);
%xlswrite(filenamelog, Matrxlog); %Modified by JV 18/12 2015, excel file is easier to work with
% Combined = [[Param; CI] AmountofMolecules];
% xlswrite(filenamerestdata,Combined);
%end

%save D coefficient data in a .txt-file
%open the output file for writing
%[Dnumber,Dxout] = hist(D,rangeD);
% fid = fopen(['Dcoeff' '_thresh' num2str(params.localizationThresh) '_win' num2str(params.trackParams.maxDisp) '_mem' num2str(params.trackParams.mem) '_nMol' num2str(nMolecules) '_MinSteps' num2str(DhistMinSteps) '_nTracks' num2str(length(D(:,1))) '.txt'],'w');
% for ina1=1:length(D(:,1))
%   fprintf(fid,'%6.3f\n',D(ina1,1)); %D coefficients  
%   %fprintf(fid,'%6.3f\n',Dxout(1,ina1)); %D range 
% end
% fclose(fid);
end
end


