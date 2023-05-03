function [Dfit,Dfitonestep,Dfittwostep,startD,tracks,pos] = SimulationLocalizationandConfinement(input, plothist)
% Simulates the presence of maximally two different particles with each
% maximally two different binding modes and a permanently bound fraction
Nparticles = input.Nparticles; % Number of particles to be simulated
Steptime = input.Steptime; % Simulates the movement per Steptime, and only checks once every frametime where the molecule is
Frametime = input.frametime;
%Steptime = Steptime*Frametime/0.01; 
%profile on
koff1_A = input.koff1_A; 
koff2_A = input.koff2_A; 
kon1_A = input.kon1_A; 
kon2_A = input.kon2_A; 
Dfree_A = input.Dfree_A;

fractionB = input.fractionB;
koff1_B = input.koff1_B; 
koff2_B = input.koff2_B; 
kon1_B = input.kon1_B; 
kon2_B = input.kon2_B; 
Dfree_B = input.Dfree_B;

sigmaerror = input.sigmaerror;
sigmasigma = 0;
%immobilefraction = input.immobilefraction;% Set what percentage of signal is immobile

koff1_A = koff1_A*Steptime;
kon1_A = kon1_A*Steptime;
koff2_A = koff2_A*Steptime;
kon2_A = kon2_A*Steptime;

koff1_B = koff1_B*Steptime;
kon1_B = kon1_B*Steptime;
koff2_B = koff2_B*Steptime;
kon2_B = kon2_B*Steptime;
confinement = input.confinement;
NFrames = input.NumberofFrames;   % Number of frames
Dim_A1 = 1e-19; %has to be different value from Dim_A1, therefore very low
Dim_A2 = 1e-20; % has to be different value from Dim_A2, therefore very low
Dim_B1 = 1e-21;%has to be different value from Dim_B1, therefore very low
Dim_B2 = 1e-22;% has to be different value from Dim_B2, therefore very low
Dfree_A = sqrt(2*Dfree_A*Steptime);
Dfree_B = sqrt(2*Dfree_B*Steptime);


%Fractions of each state depending on their equilibrium values from koff and kon
c2_A = 1/(1+(kon1_A/koff1_A)+(kon2_A/koff2_A)); 
c1_A = (kon1_A/koff1_A)*c2_A;

c2_B = 1/(1+(kon1_B/koff1_B)+(kon2_B/koff2_B));
c1_B = (kon1_B/koff1_B)*c2_B;



Numpop = 3;
rng('shuffle', 'simdTwister');

%Tracking parameters
trackparams.maxDisp = input.trackingwindow;% Maximum tracking window
trackparams.maxDisp = trackparams.maxDisp*input.pixelsize;
trackparams.mem = 0;      % Memory
trackparams.dim = 2;
trackparams.good = 0;
trackparams.quiet = 1;
fulltrack = 0; % if 1 uses the 'standard tracking software' which is far slower but can take things like shot noise and density into account
fullhistogram = 0;% if 1 uses the 'standard histogram software' which is far slower but can fit different states and take shot noise etc into account

%Interchanging parameters
interchangingstate = true; % if you want to have states interchanging, set to true
recordstates = 0; % Parameter controls whether you want to record the states over time (which you can use later for (1 active))

plotmix = 0; % If you want to plot the simulated populations 
threshold = 0.9; % Parameter controls what percentage of the time it should be in this state to be assigned to the state in the histogram

%Size parameters and additional simulations
xdimcell = input.lengthcell; %length of cell in um 
radiusofcell = input.radiusofcell; %um
runtime = Frametime/10; %s
densitylocalizations = 0; % if you want to introduce bias by having too many fluorophores activated at same time (if 0 not active)
densityshotnoise = 0;% if you want to introduce effects of shot noise on data(if 0 not active)

% Autofluorescence parameters 
autofluorescence = 1;% If you want to include immobile autofluorescent particles set to 1

%% Part to simulate movement in cell and confinement effects 
l = 1; % position index 
t = 0; %start time
cont = true(Nparticles,1);
Poslist = zeros(Nparticles*(NFrames+1),3);
startpos=zeros(Nparticles,3);
sumcont = sum(cont(:));
%% Give each particle a starting position
if confinement == 1
    while sumcont > 0

        startpos(cont,1)=xdimcell*rand(sumcont, 1,'single');
      startpos(cont,2:3,1) = (rand(sumcont,2,'single')-0.5)*(radiusofcell*2);
      squared = startpos(cont,2).^2 + startpos(cont,3).^2;
      capcond1 = ...
             ...
            ((startpos(cont,1) < radiusofcell) .*...
            (startpos(cont,1)-radiusofcell).^2 +squared <radiusofcell^2);
    capcond2 = ...
             ...
            ((startpos(cont,1) >= radiusofcell) .*...    
            (startpos(cont,1)>(xdimcell-radiusofcell)) .* ...
             (startpos(cont,1)-(xdimcell-radiusofcell)).^2 + squared <radiusofcell^2);
    capcond = capcond1 & capcond2;
    particlesincell(cont,1) = startpos(cont,1)>0 & (startpos(cont,1)<xdimcell) & (squared<radiusofcell^2) & capcond;%.*running_particles;
    cont(cont,1) = cont(cont,1)-particlesincell(cont,1);
    sumcont = sum(cont);
    end
else
     startpos(cont,1:3)=rand(sumcont, 3,'single');
end
%% Give each particle a starting Diffusion coefficient value
startD = rand(Nparticles,1);
D(1:Nparticles,1)= Dim_A2;
D(startD<(c1_A+c2_A))= Dfree_A;
D(startD<(c1_A))= Dim_A1;
startD = D;

% Subunit only fractions from previously calculated 
if fractionB > 0
startDCas8e = rand(Nparticles,1);
Cas8eparticle = rand(Nparticles,1);
Cas8e_3 = (startDCas8e>(c1_B+c2_B)).*(Cas8eparticle<fractionB);
Cas8e_2 = (startDCas8e<(c1_B+c2_B)).*(Cas8eparticle<fractionB);
Cas8e_1 = (startDCas8e<c1_B).*(Cas8eparticle<fractionB);
D(Cas8e_3==1,1)= Dim_B2;
D(Cas8e_2==1,1)= Dfree_B;
D(Cas8e_1==1,1)= Dim_B1;
end

startD = D;
ll = 1;
while t < (runtime + (NFrames+0.6)*Frametime)    
%% Interchanging part
if interchangingstate == true

change = rand(Nparticles,1);

% Change free to bound and bound to free depending on probability
FreetoBound1 = (D==Dfree_A)&(change<kon1_A);
FreetoBound2 = (D==Dfree_A)&(change>(1-kon2_A));
BoundtoFree1 = (D ==Dim_A1)&(change<koff1_A);
BoundtoFree2 = (D ==Dim_A2)&(change<koff2_A);

% Change diffusion coefficients of those particles that changed states
D(FreetoBound1) = Dim_A1;
D(FreetoBound2) = Dim_A2;
D(BoundtoFree1) = Dfree_A;
D(BoundtoFree2) = Dfree_A;

if fractionB > 0
FreetoBound1Cas8e = (D==Dfree_B)&(change<kon1_B);
FreetoBound2Cas8e = (D==Dfree_B)&(change>(1-kon2_B));
BoundtoFree1Cas8e = (D ==Dim_B1)&(change<koff1_B);
BoundtoFree2Cas8e = (D ==Dim_B2)&(change<koff2_B);
D(FreetoBound1Cas8e) = Dim_B1;
D(FreetoBound2Cas8e) = Dim_B2;
D(BoundtoFree1Cas8e) = Dfree_B;
D(BoundtoFree2Cas8e) = Dfree_B;
end

% Part to record states so you can later follow the states and determine
% fraction of time spent in each state for each particle
    if recordstates == 1
        startD = [startD D];
    end
end

% make stepsize based on D
%k = sqrt(2*D*Steptime);
k = D;
% initialise loop
running_particles = true(Nparticles,1);
%% Set movement
%mostly based on Koen Martens script, however here the particles that are not 
%inside of the cell after step are first all resolved before moving to next time point

while sum(running_particles>0)
%         %k(running_particles == 1) = sqrt(2.*D(running_particles==1)*Steptime);
% % Change in xyz position
% dxyz(running_particles ==1,:) = k(running_particles == 1).*(randn(sum(running_particles),3));
% %get the new position 
% %newpos(running_particles ==1,:) = startpos(running_particles ==1,:)+dxyz(running_particles ==1,:);
% newpos = running_particles.*startpos+dxyz;
% %Check if the molecules are inside the cell
% squared = newpos(:,2).^2 + newpos(:,3).^2;
% % capcond == 1 is necessary for continuation
% %Splitting up in two parts 
% capcond1 = ...
%         running_particles .* ...
%         ((newpos(:,1) < radiusofcell) .*...
%         (newpos(:,1)-radiusofcell).^2 +squared <radiusofcell^2);
% capcond2 = ...
%         running_particles .* ...
%         ((newpos(:,1) >= radiusofcell) .*...    
%         (newpos(:,1)>(xdimcell-radiusofcell)) .* ...
%          (newpos(:,1)-(xdimcell-radiusofcell)).^2 + squared <radiusofcell^2);
% capcond = capcond1 .* capcond2;
%      
% %now check all particles that are in the cell at the end of the step
% particlesincell = (newpos(:,1)>0) .* (newpos(:,1)<xdimcell) .* (squared<radiusofcell^2) .* capcond;%.*running_particles;
% startpos(particlesincell==1,:) = newpos(particlesincell==1,:);
% %startpos = newpos.*particlesincell + startpos.*(~particlesincell); 
% running_particles = running_particles - particlesincell;
% %     k(running_particles == 1) = sqrt(2.*D(running_particles==1)*Steptime);
sumrun = sum(running_particles); 
% Change in xyz position
if sumrun == Nparticles
dxyz = k.*(randn(sumrun,3));
newpos = startpos+dxyz;
else
dxyz = k(running_particles).*(randn(sumrun,3));
newpos = startpos(running_particles,:)+dxyz;  
end
%get the new position 
%newpos(running_particles ==1,:) = startpos(running_particles ==1,:)+dxyz(running_particles ==1,:);
%newpos = [];
if confinement == 1
%Check if the molecules are inside the cell
squared = newpos(:,2).^2 + newpos(:,3).^2;
% capcond == 1 is necessary for continuation
%Splitting up in two parts 
capcond1 = ...
         ...
        ((newpos(:,1) < radiusofcell) .*...
        (newpos(:,1)-radiusofcell).^2 +squared <radiusofcell^2);
capcond2 = ...
         ...
        ((newpos(:,1) >= radiusofcell) .*...    
        (newpos(:,1)>(xdimcell-radiusofcell)) .* ...
         (newpos(:,1)-(xdimcell-radiusofcell)).^2 + squared <radiusofcell^2);
capcond = capcond1 & capcond2;     
%now check all particles that are in the cell at the end of the step
particlesincell = newpos(:,1)>0 & (newpos(:,1)<xdimcell) & (squared<radiusofcell^2) & capcond;%.*running_particles;
else
particlesincell = ones(numel(running_particles),1);
end
if sumrun == Nparticles
newpos(~particlesincell,:) = startpos(~particlesincell,:);
startpos = newpos;
%    startpos(particlesincell,:) = newpos(particlesincell,:);
running_particles = running_particles & ~particlesincell;
else
particlesincell2 = false(Nparticles,1);
particlesincell2(running_particles)=particlesincell;
startpos(particlesincell2,:) = newpos(particlesincell,:);
running_particles = running_particles & ~particlesincell2;
end

%startpos = newpos.*particlesincell + startpos.*(~particlesincell); 


if sum(running_particles)>0
%newpos = [];
newpos2 = newpos(~particlesincell,:)-2*dxyz(~particlesincell,:);
%Check if the molecules are inside the cell
squared2 = newpos2(:,2).^2 + newpos2(:,3).^2;
% capcond == 1 is necessary for continuation
%Splitting up in two parts 
capcond12 = ...
         ...
        ((newpos2(:,1) < radiusofcell) .*...
        (newpos2(:,1)-radiusofcell).^2 +squared2 <radiusofcell^2);
capcond22 = ...
         ...
        ((newpos2(:,1) >= radiusofcell) .*...    
        (newpos2(:,1)>(xdimcell-radiusofcell)) .* ...
         (newpos2(:,1)-(xdimcell-radiusofcell)).^2 + squared2 <radiusofcell^2);
capcond3 = capcond12 & capcond22;     
%now check all particles that are in the cell at the end of the step
particlesincell = newpos2(:,1)>0 & (newpos2(:,1)<xdimcell) & (squared2<radiusofcell^2) & capcond3;%.*running_particles;
particlesincell2 = false(Nparticles,1);
particlesincell2(running_particles)=particlesincell;
startpos(particlesincell2,:) = newpos2(particlesincell,:);
%startpos = newpos.*particlesincell + startpos.*(~particlesincell); 
running_particles = running_particles & ~particlesincell2;    
end

end

%% Displace + possibly record location
% if newpos(1)>0 && newpos(1)<xdimcell && squared < radiusofcell^2 && capcond == 1
%     startpos = newpos;
try
    t = t+Steptime;
    if  t > (runtime-0.5*Frametime) && (mod(t-runtime,Frametime) < 0.0000000002 || mod(t-runtime,Frametime) > Frametime-0.0000000002) 
    %if  t > (runtime-0.5*Frametime) && (mod(t-runtime,Frametime) < 0.00000002 || mod(t-runtime,Frametime) < 0.00000002) 
    Poslist(l:(NFrames+1):end)= startpos;
    l = l+1; 
    end
catch
    keyboard
end
end
%% Add localization error + add recorded times depending on density of localizations
Possize = size(Poslist);
% Can add localizationerror but also can have a certain distribution of
% brightness (localization error) for each particle (sigmasigma). 
sigmaerror = abs(normrnd(sigmaerror,sigmasigma,Nparticles,1));
sigmaerror = repelem(sigmaerror,NFrames+1);
sigmaerror = repmat(sigmaerror,1,3);
localizationerror = normrnd(0,sigmaerror);
Poslist = Poslist+ localizationerror;
k = 0;
% Adds time points in case you want to investigate the effect of density of
% localizations

for i = 1:numel(Poslist(:,1))
Poslist(i,4)=i+k;
if mod(i,NFrames+1) == 0  
if densitylocalizations > 0   
    nexttime = round(-log(rand)/densitylocalizations);
    k = k +nexttime;
    else
    k = k+2;
    end
end
%else
%diff = setdiff(1:2*numel(Poslist(:,1)),NFrames+2:NFrames+3:2*numel(Poslist(:,1))); 
%diff = setdiff(diff,NFrames+3:NFrames+3:2*numel(Poslist(:,1)));
%Poslist(:,4) = diff(1:numel(Poslist(:,1)));
end

%% Can add shot noise if needed
if densityshotnoise > 0
ll = 0;
while ll < max(Poslist(:,4))
nexttime = round(-log(rand)/densityshotnoise);
ll = ll+nexttime;
while cont == 0
startpos(1) = xdimcell*rand(1);
startpos(2:3,1) = (rand(2,1)-radiusofcell)/(radiusofcell*2);
if startpos(2)^2+startpos(3)^2<radiusofcell^2
    if startpos(1)<radiusofcell
        if(startpos(1)-radiusofcell)^2 + startpos(2)^2+startpos(3)^2<radiusofcell^2
        cont = 1;
        end
    elseif startpos(1)>(xdimcell-radiusofcell)
        if(startpos(1)-(xdimcell-radiusofcell))^2 + startpos(2)^2+startpos(3)^2<radiusofcell^2
        cont = 1;
        end
    else
    cont = 1;
    end
end
end
Poslist = [Poslist; [startpos(1) startpos(2) startpos(3) ll]];
end
end  
%% Prepare for tracking and histograms
% Initialisations of position list
pos(:,1:2) = Poslist(:,1:2);
% pos(:,1:2) = pos(:,1:2)-min(pos(:))+1;
pos(:,3) = Poslist(:,4);
MSD = zeros(Nparticles*(NFrames+1),1);
if fulltrack == 1 % the full tracking script
    tracks = trackWithDummy(pos,trackparams);
else 
    % The fast shortcut script for tracking which can still keep tracking
    % window into account
    MSD(2:end) = sqrt(diff(pos(:,1)).^2+diff(pos(:,2)).^2);
    window = MSD>trackparams.maxDisp;
    window(NFrames+2:NFrames+1:end)=1;
    track = 1 + cumsum(window);
    tracks = [pos track];
%     track(1) = 1;
%     l = 1;
%     m = 1;
%     for i = 2:numel(pos(:,1))
%         if mod(i-1,NFrames+1) == 0
%             m =1;
%             l = l+1;
%         else
%             MSD(i) = sqrt((pos(i,1)-pos(i-1,1))^2+(pos(i,2)-pos(i-1,2))^2);
%             if MSD(i) >trackparams.maxDisp % disconnects localizations that are further apart then tracking window
%                 l = l+1;
%             end
%         end
%         track(i)=l;
%     end
%     track = track';
%     tracks = [pos track];
end
% This part gives the diffusion coefficients, can deal with both tracking
% versions (fast and slow). For plotting and fitting and also if you introduce shot noise data, use fullhistogram
% = 1, if you only need the Dfit distribution then you can use the faster
% script below
if fullhistogram == 1
        Dfit = histDadapted(tracks, Frametime, Numpop, plothist);
else    
MSD = MSD.^2;
MSD = reshape(MSD,[(NFrames+1),Nparticles])';
MSD = MSD(:,2:end);
MSD(MSD>trackparams.maxDisp^2)=NaN;


MSDav = mean(MSD,2);
MSDav = MSDav(~isnan(MSDav));
MSD(isnan(MSDav),:) = [];
Dfitonestep = MSD(:,1:2:end)/(4*Frametime);
Dfitonestep = Dfitonestep(:)';

if NFrames==4
Dfittwostep =  mean(MSD(:,1:2),2);
%Dfittwostep = [Dfittwostep; mean(MSD(:,3:4),2)];
Dfittwostep = Dfittwostep'/(4*Frametime);
else
Dfittwostep = 3;
end
sqrtMSDav = sqrt(MSDav);
Dfit = MSDav/(4*Frametime);
%Dfit = Dfit(Dfit>0);
Dfit = Dfit';

end
%% Part to plot the different populations in a histogram

if plotmix == 1
% only pick the part of the states that is actually recorded
startD = startD(:,runtime/Steptime:(runtime/Steptime+(NFrames*Frametime/Steptime)));
% Find how often each state occurs in the total dataset
[rowsD1, ~] = find(startD==Dim_A1);
[rowsD2, ~] = find(startD==Dfree_A);
[rowsD3, ~] = find(startD==Dim_A2);

% Find how often each state occurs in the total dataset
occurD1 = histc(rowsD1,1:Nparticles);
occurD2 = histc(rowsD2,1:Nparticles);
occurD3 = histc(rowsD3,1:Nparticles);
columns= size(startD);
columns = columns(2);

% Find the diffusion coefficients of each state
% The threshold value determines what fraction of the time it needs to have
% the diffusion coefficient of that state
State1 = occurD1>threshold*columns;
State2 = occurD2>threshold*columns;
State3 = occurD3>threshold*columns;
mix = ~(State1+State2+State3);
MSDstate1 = MSDav(State1);
MSDstate2 = MSDav(State2);
MSDstate3 = MSDav(State3);
MSDstatemix = MSDav(mix);
Dstate1 =  MSDstate1(MSDstate1>0)/(4*Frametime);
Dstate2 =  MSDstate2(MSDstate2>0)/(4*Frametime);
Dstate3 =  MSDstate3(MSDstate3>0)/(4*Frametime);
Dstatemix =  MSDstatemix(MSDstatemix>0)/(4*Frametime);

% Plot the different states as histogram
figure
hold on
axis([0 8 0 3000])
histogram(Dstate1,[0:0.1:8])
histogram(Dstate2,[0:0.1:8])
histogram(Dstatemix,[0:0.1:8])

figure

% Plot the different states as lines
figure
hold on
%axis([0 8 0 3000])
counts1 = histc(Dstate1,[0:0.1:8]);
plot([0:0.1:8],10*counts1/sum(Dfit));
counts2 = histc(Dstate2,[0:0.1:8]);
plot([0:0.1:8],10*counts2/sum(Dfit));
counts3 = histc(Dstatemix,[0:0.1:8]);
plot([0:0.1:8],10*counts3/sum(Dfit));
figure
counts = [counts1' counts2' counts3'];
bar(counts,'stacked')
end
