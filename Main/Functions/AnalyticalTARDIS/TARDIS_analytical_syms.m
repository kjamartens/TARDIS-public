function [TARDIS_syms] = TARDIS_analytical_syms(settings)
%% INFO
%First run GenerateTARDISAnalyticalSettings, then this.

%Obtain relative lambda_bleach
lambda_bleach_movLength = settings.lambda_bleach/settings.movLength;
%Prevent Inf-artefacts
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
    figure(2);clf(2);
end
%% Inter-particle distance distribution
syms x

eq1 = 2*x*(-4*sqrt(x^2)+pi()+x^2);
eq2 = 2*x*(-2+4*asin(1/sqrt(x^2))+4*sqrt(x^2-1)-pi-x^2);
InterParticleDistDistributionSyms = ...
    piecewise(x<=1,eq1,x>1,eq2);

if settings.figureCreation
    subplot(3,6,[1,2,3])
    hold on
    title('Inter-particle distance distribution')    
    fplot(InterParticleDistDistributionSyms,[0 1.41],'k-','LineWidth',2,'DisplayName','Probability density')
    xlabel('Distance between two random points on a 1-by-1 plane (AU)')
    ylabel('pdf')
end

clear x
%% Intra-particle distance distribution
syms dt tau D sigma st dt x dim
%Formalise sub-parts of the formula
dttauSyms = dt*tau;
DandsigSyms = (D+sigma^2/dttauSyms);
RSyms = st/dt*(1/6);
mainFactorSyms = exp(-x^2/(2*(dim*DandsigSyms*dttauSyms-2*RSyms*DandsigSyms*dttauSyms)));
normFactorSyms = x/((dim*DandsigSyms*dttauSyms-2*RSyms*DandsigSyms*dttauSyms));

pdf_fSyms = mainFactorSyms.*normFactorSyms;

for tauv = settings.tauRange_distribution
    tau = tauv;
    D = settings.D;
    sigma = settings.sigma;
    st = settings.st;
    dt = settings.dt;
    dim = settings.dim;
    ndata = 1/settings.dx;
    pdf_f_subs(tau) = subs(pdf_fSyms);
    pdf_f_subs_Norm(tau) = subs(pdf_fSyms/ndata);
end
if settings.figureCreation == 1
    subplot(3,6,[4,5,6])
    hold on
    title('Intra-particle distance distribution')
    for tau = settings.tauRange_distribution
        fplot(pdf_f_subs_Norm(tau),[0 0.2],'DisplayName',['\tau = ' num2str(tau) ' frames'])
    end
    axis([ 0 0.2 0 inf])
    legend()
    xlabel('Distance between two \tau-spaced trajectory points on a 1-by-1 plane (AU)')
    ylabel('Probability')
end
clear dt tau D sigma st dt x dim
%% Inter-particle contribution
syms x d_spurious d_traj movLength lambda_bleach p_nonLocalized lambda_bright lambda_dark tau

P_bright = lambda_bright/(lambda_bright+lambda_dark);

n_locs_f_traj = ((2^(1/lambda_bleach))/((2^(1/lambda_bleach)-1)))*d_traj*P_bright;
n_locs_f_trajSyms = n_locs_f_traj*(1-p_nonLocalized);
n_locs_spSyms = (d_spurious)*(1-p_nonLocalized);
n_locs_f_totSyms = n_locs_f_trajSyms+n_locs_spSyms;

if settings.figureCreation
    d_traj = settings.d_traj;
    movLength = settings.movLength;
    lambda_bleach = settings.lambda_bleach;
    p_nonLocalized = settings.p_nonLocalized;
    d_spurious = settings.d_spurious;
    lambda_bright = settings.lambda_bright;
    lambda_dark = settings.lambda_dark;

    subplot(3,6,[7,8,9])
    hold on
    title('Inter-particle contribution')
    fplot(subs(n_locs_f_trajSyms),[1 movLength],'r-','DisplayName','Trajectory-based')
    fplot(subs(n_locs_spSyms),[1 movLength],'b-','DisplayName','Spurious-based')
    fplot(subs(n_locs_f_totSyms),[1 movLength],'k-','DisplayName','Total')
    xlabel('Frame')
    ylabel('Number of localizations')
    legend()
    axis([1, movLength, 0, eval(subs(n_locs_f_totSyms))*1.1])
end

clear d_spurious d_traj movLength lambda_bleach n_tracks p_nonLocalized
%% Intra-particle contribution
syms tau lambda_bleach p_nonLocalized x lambda_bright lambda_dark

n_intraParticle_p1 = ((2.^((1-tau)./lambda_bleach))/(2.^(1/lambda_bleach)-1));
n_intraParticle_p2 = ((lambda_bright.*(lambda_dark.*2.^(-tau./lambda_dark-tau./lambda_bright)+lambda_bright))./((lambda_dark+lambda_bright)^2));
n_intraParticle_p3 = (1-p_nonLocalized)^2;
n_intraParticleSyms = n_intraParticle_p1.*n_intraParticle_p2.*n_intraParticle_p3;

if settings.figureCreation
    lambda_bleach = settings.lambda_bleach;
    p_nonLocalized = settings.p_nonLocalized;
    lambda_bright = settings.lambda_bright;
    lambda_dark = settings.lambda_dark;
    tau = x; %We plot over tau
    subplot(3,6,[10,11,12])
    hold on
    title('Intra-particle contribution')
    fplot(subs(n_intraParticleSyms),[1 max(settings.tauRange_contribution)],'k-','DisplayName','n_{intraParticle}')
    xlabel('Temporal shift (\tau) (frames)')
    ylabel('Number of links per track')
    legend()
    axis([1, max(settings.tauRange_contribution), 0, max(eval(subs(subs(n_intraParticleSyms),tau,1)))*1.1])
end

%% Creation of 'full TARDIS result'
syms dt tau D sigma st dt x dim d_spurious d_traj movLength lambda_bleach p_nonLocalized lambda_bright lambda_dark

%Inter-particle distance distribution
eq1 = 2*x*(-4*sqrt(x^2)+pi()+x^2);
eq2 = 2*x*(-2+4*asin(1/sqrt(x^2))+4*sqrt(x^2-1)-pi-x^2);
InterParticleDistDistributionSyms = ...
    piecewise(x<=1,eq1,x>1,eq2);

%Intra-particle distance distribution
%Formalise sub-parts of the formula
dttauSyms = dt*tau;
DandsigSyms = (D+sigma^2/dttauSyms);
RSyms = st/dt*(1/6);
mainFactorSyms = exp(-x^2/(2*(dim*DandsigSyms*dttauSyms-2*RSyms*DandsigSyms*dttauSyms)));
normFactorSyms = x/((dim*DandsigSyms*dttauSyms-2*RSyms*DandsigSyms*dttauSyms));
IntraParticleDistDistributionSyms = mainFactorSyms.*normFactorSyms;

%Inter-particle contribution
P_bright = lambda_bright/(lambda_bright+lambda_dark);

n_locs_f_traj = ((2^(1/lambda_bleach))/((2^(1/lambda_bleach)-1)))*d_traj*P_bright;
n_locs_f_trajSyms = n_locs_f_traj*(1-p_nonLocalized);
n_locs_spSyms = (d_spurious)*(1-p_nonLocalized);
n_locs_f_totSyms = n_locs_f_trajSyms+n_locs_spSyms;

%Intra-particle contribution 
n_intraParticle_p1 = ((2.^((1-tau)./lambda_bleach))/(2.^(1/lambda_bleach)-1));
n_intraParticle_p2 = ((lambda_bright.*(lambda_dark.*2.^(-tau./lambda_dark-tau./lambda_bright)+lambda_bright))./((lambda_dark+lambda_bright)^2));
n_intraParticle_p3 = (1-p_nonLocalized)^2;
IntraParticleContributionSyms = n_intraParticle_p1.*n_intraParticle_p2.*n_intraParticle_p3;

%Calculate the fractions of inter/intra-particles at all taus
n_tracksSyms = d_traj*movLength;
fraction_interParticleSyms = ((movLength-tau)*n_locs_f_totSyms*n_locs_f_totSyms-IntraParticleContributionSyms*n_tracksSyms)/((movLength-tau)*n_locs_f_totSyms*n_locs_f_totSyms);
fraction_intraParticleSyms = 1-fraction_interParticleSyms;

TotalInterParticleDistributionSyms = fraction_interParticleSyms*InterParticleDistDistributionSyms;
TotalIntraParticleDistributionSyms = fraction_intraParticleSyms*IntraParticleDistDistributionSyms;
TotalTARDISDistributionSyms = TotalInterParticleDistributionSyms+TotalIntraParticleDistributionSyms;

%Store for export from this function
TARDIS_syms = TotalTARDISDistributionSyms;

for tauv = settings.tauRange_distribution
    if tauv <=3
        if settings.figureCreation
            subplot(3,6,[13+(tauv-1)*2, 14+(tauv-1)*2])

            dt = settings.dt;
            tau = tauv;
            D = settings.D;
            sigma = settings.sigma;
            st = settings.st;
            dt = settings.dt;
            dim = settings.dim;
            d_spurious = settings.d_spurious;
            d_traj = settings.d_traj;
            movLength = settings.movLength;
            lambda_bleach = settings.lambda_bleach;
            p_nonLocalized = settings.p_nonLocalized;
            lambda_bright = settings.lambda_bright;
            lambda_dark = settings.lambda_dark;

            hold on
            fplot(subs(TotalInterParticleDistributionSyms),[0 max(settings.xp)],'b-','LineWidth',2,'DisplayName','Inter-particle')
            fplot(subs(TotalIntraParticleDistributionSyms),[0 max(settings.xp)],'r-','LineWidth',2,'DisplayName','Intra-particle')
            pl{tauv} = fplot(subs(TotalTARDISDistributionSyms),[0 max(settings.xp)],'k-','LineWidth',2,'DisplayName','Combination');
            legend()
            title(['TARDIS distribution at \tau = ' num2str(tauv) ' frames'])
            xlabel('Jump distance (AU)')
            ylabel('Probability')
        end
    end

end

if settings.figureCreation
    %Normalise axis between all three
    maxyval = 0;
    for tau = 1:3
        if max(pl{tau}.YData) > maxyval
            maxyval = max(pl{tau}.YData);
        end
    end
    for tau = 1:3
        subplot(3,6,[13+(tau-1)*2, 14+(tau-1)*2])
        axis([-inf, inf, 0, 1.1*maxyval])
    end
end

clear dt tau D sigma st dt x dim d_spurious d_traj movLength lambda_bleach p_nonLocalized
end