function [distribution_TARDIS,distribution_interParticle,distribution_intraParticle,fraction_interParticle,xdata,n_linkages] = TARDIS_analytical(settings)
%% INFO
%First run GenerateTARDISAnalyticalSettings, then this.

%Obtain relative lambda
lambda_bleach_movLength = settings.lambda_bleach/settings.movLength;

if settings.lambda_bright == Inf
    settings.lambda_bright = 1e99;
end
if settings.lambda_dark == Inf
    settings.lambda_dark = 1e99;
end
%Prevent division-by-zero artefacts
if settings.lambda_bright == 0
    settings.lambda_bright = 1e-99;
end
if settings.lambda_dark == 0
    settings.lambda_dark = 1e-99;
end

%% Figure initialisation
if settings.figureCreation
    figure(1);clf(1);
end
%% Inter-particle distance distribution
dd1 = [0:settings.dx:1];
dd2 = [1:settings.dx:sqrt(2)];

if settings.figureCreation
    subplot(3,6,[1,2,3])
    hold on
    title('Inter-particle distance distribution')
    plot(dd1,2.*dd1.*(-4*sqrt(dd1.^2)+pi()+dd1.^2),'k-','LineWidth',2,'DisplayName','Probability density')
    plot(dd2,2.*dd2.*(-2+4.*asin(1./sqrt(dd2.^2))+4.*sqrt(dd2.^2-1)-pi()-dd2.^2),'k-','LineWidth',2,'DisplayName','Probability density')
    xlabel('Distance between two random points on a 1-by-1 plane (AU)')
    ylabel('pdf')
end
%% Intra-particle distance distribution
pdf_f = zeros(size(settings.tauRange_distribution,1),size(settings.xp,2));
for tau = settings.tauRange_distribution
    %Formalise sub-parts of the formula
    dttau = settings.dt*tau;
    Dandsig = (settings.D+settings.sigma^2/dttau);
    R = settings.st/settings.dt*(1/6);
    %Create the formula
    mainfactor = exp(-settings.xp.^2./(2*(settings.dim*Dandsig*dttau-2*R*Dandsig*dttau))); %Main factor
    normfactor = settings.dx.*settings.xp ./ ((settings.dim*Dandsig*dttau-2*R*Dandsig*dttau));
    %Normalize
    pdf_f(tau,:) = mainfactor.*normfactor; %Calculate
end
if settings.figureCreation == 1
    subplot(3,6,[4,5,6])
    hold on
    title('Intra-particle distance distribution')
    for tau = settings.tauRange_distribution
        plot(settings.xp,pdf_f(tau,:),'DisplayName',['\tau = ' num2str(tau) ' frames'])
    end
    axis([ 0 0.2 0 inf])
    legend()
    xlabel('Distance between two \tau-spaced trajectory points on a 1-by-1 plane (AU)')
    ylabel('Probability')
end

%% Inter-particle contribution
P_open_form = settings.lambda_bright/(settings.lambda_bright+settings.lambda_dark);

n_tracks = settings.d_traj*settings.movLength;
n_locs_f_traj = ((2^(1/settings.lambda_bleach))/((2^(1/settings.lambda_bleach)-1)))*settings.d_traj*P_open_form;
n_locs_f_traj = n_locs_f_traj*(1-settings.p_nonLocalized);
n_locs_sp = (settings.d_spurious)*(1-settings.p_nonLocalized);
n_locs_f_tot = n_locs_f_traj+n_locs_sp;

if settings.figureCreation
    subplot(3,6,[7,8,9])
    hold on
    title('Inter-particle contribution')
    plot([1,settings.movLength],[n_locs_f_traj, n_locs_f_traj],'r-','DisplayName','Trajectory-based')
    plot([1,settings.movLength],[n_locs_sp,n_locs_sp],'b-','DisplayName','Spurious-based')
    plot([1,settings.movLength],[n_locs_f_tot,n_locs_f_tot],'k-','DisplayName','Total')
    xlabel('Frame')
    ylabel('Number of localizations')
    legend()
    axis([1,settings.movLength, 0, max(n_locs_f_tot)*1.1])
end

%% Intra-particle contribution
%Partials
n_intraParticle_p1 = ((2.^((1-settings.tauRange_contribution)./settings.lambda_bleach))/(2.^(1/settings.lambda_bleach)-1));
n_intraParticle_p2 = ((settings.lambda_bright.*(settings.lambda_dark.*2.^(-settings.tauRange_contribution./settings.lambda_dark-settings.tauRange_contribution./settings.lambda_bright)+settings.lambda_bright))./((settings.lambda_dark+settings.lambda_bright)^2));
n_intraParticle_p3 = (1-settings.p_nonLocalized)^2;
%Full
n_intraParticle = n_intraParticle_p1.*n_intraParticle_p2.*n_intraParticle_p3;

if settings.figureCreation
    subplot(3,6,[10,11,12])
    hold on
    title('Intra-particle contribution')
    plot(settings.tauRange_contribution,n_intraParticle,'k-','DisplayName','n_{intraParticle}')
    xlabel('Temporal shift (\tau) (frames)')
    ylabel('Number of links per track')
    legend()
    axis([1, max(settings.tauRange_contribution), 0, max(n_intraParticle)*1.1])
end
%% Creation of 'full TARDIS result'
%Initiate arrays
distribution_interParticle = zeros(size(settings.xp,2),size(settings.tauRange_distribution,2));
distribution_intraParticle = zeros(size(settings.xp,2),size(settings.tauRange_distribution,2));
distribution_TARDIS = zeros(size(settings.xp,2),size(settings.tauRange_distribution,2));
fraction_interParticle = zeros(0,size(settings.tauRange_distribution,2));
fraction_intraParticle = zeros(0,size(settings.tauRange_distribution,2));
n_linkages = zeros(size(settings.tauRange_distribution,2),2);

%Loop over wanted tau ranges
for tau = settings.tauRange_distribution
    if tau <=3
        if settings.figureCreation
            subplot(3,6,[13+(tau-1)*2, 14+(tau-1)*2])
        end
    end

    %Number of jumps between f and f+tau FOR A SINGLE TRAJECTORY
    n_intraParticle_p1 = ((2.^((1-tau)./settings.lambda_bleach))/(2.^(1/settings.lambda_bleach)-1));
    n_intraParticle_p2 = ((settings.lambda_bright.*(settings.lambda_dark.*2.^(-tau./settings.lambda_dark-tau./settings.lambda_bright)+settings.lambda_bright))./((settings.lambda_dark+settings.lambda_bright)^2));
    n_intraParticle_p3 = (1-settings.p_nonLocalized)^2;
    n_intraParticle = n_intraParticle_p1.*n_intraParticle_p2.*n_intraParticle_p3;

    %Inter-particle
    P_open_form = settings.lambda_bright/(settings.lambda_bright+settings.lambda_dark);
    
    n_tracks = settings.d_traj*settings.movLength;
    n_locs_f_traj = ((2^(1/settings.lambda_bleach))/((2^(1/settings.lambda_bleach)-1)))*settings.d_traj*P_open_form;
    n_locs_f_traj = n_locs_f_traj*(1-settings.p_nonLocalized);
    n_locs_sp = (settings.d_spurious)*(1-settings.p_nonLocalized);
    n_locs_f_tot = n_locs_f_traj+n_locs_sp;
    
    %Calculate the fractions of inter/intra-particles
    fraction_interParticle(tau) = ((settings.movLength-tau)*n_locs_f_tot*n_locs_f_tot-n_intraParticle*n_tracks)/((settings.movLength-tau)*n_locs_f_tot*n_locs_f_tot);
    fraction_intraParticle(tau) = 1-fraction_interParticle(tau);

    %Also save the number of linkages
    n_linkages(tau,:) = [(settings.movLength-tau)*n_locs_f_tot*n_locs_f_tot,n_intraParticle*n_tracks];

    %Calculate inter-particle distribution at this x - different below 1
    %and above 1
    Xbelow1 = settings.xp < 1;
    Xabove1 = settings.xp >= 1;
    distribution_interParticle(:,tau) = Xbelow1.*(settings.dx*2.*settings.xp.*(-4*sqrt(settings.xp.^2)+pi()+settings.xp.^2))+...
        real(Xabove1.*(settings.dx*2.*settings.xp.*(-2+4.*asin(1./sqrt(settings.xp.^2))+4.*sqrt(settings.xp.^2-1)-pi()-settings.xp.^2)));
    distribution_interParticle(1,tau) = 0;
    %Formalise sub-parts of the formula
    dttau = settings.dt*tau;
    Dandsig = (settings.D+settings.sigma^2/dttau);
    R = settings.st/settings.dt*(1/6);
    %Calculate intra-particle distribution at this x
    mainfactor = exp(-settings.xp.^2./(2*(settings.dim*Dandsig*dttau-2*R*Dandsig*dttau))); %Main factor
    normfactor = settings.dx.*settings.xp ./ ((settings.dim*Dandsig*dttau-2*R*Dandsig*dttau));
    %Normalize
    distribution_intraParticle(:,tau) = mainfactor.*normfactor;
    
    %Sum together based on fractions
    distribution_TARDIS(:,tau) = distribution_interParticle(:,tau)*fraction_interParticle(tau)+distribution_intraParticle(:,tau)*fraction_intraParticle(tau);
    if tau <=3
        if settings.figureCreation
            hold on
            plot(settings.xp,distribution_interParticle(:,tau)*fraction_interParticle(tau),'b-','LineWidth',2,'DisplayName','Inter-particle')
            plot(settings.xp,distribution_intraParticle(:,tau)*fraction_intraParticle(tau),'r-','LineWidth',2,'DisplayName','Intra-particle')
            plot(settings.xp,distribution_TARDIS(:,tau),'k-','LineWidth',2,'DisplayName','Combination')
            legend()
            title(['TARDIS distribution at \tau = ' num2str(tau) ' frames'])
            xlabel('Jump distance (AU)')
            ylabel('Probability')
        end
    end
end
if settings.figureCreation
    %Normalise axis between all three
    for tau = 1:3
        subplot(3,6,[13+(tau-1)*2, 14+(tau-1)*2])
        axis([-inf, inf, 0, 1.1*max(max(distribution_TARDIS))])
    end
end
%Only for output, create xdata
xdata = settings.xp';
end