%% Text output from TARDIS settings
%Goal: get human-readable text output from TARDIS settings for use in MM
%etc

inputloc = '\\ifmb-nas\ag endesfelder\Data\Koen\Articles\2020_relativeDisplacementTracking\SIData_figure_TrackingChallenge\TARDIS_Results\MICROTUBULE snr 1 density high_TARDISresult.mat';

load(inputloc);
settings = appinfo.settingsURDA;

if settings.createJDsonly
    out1 = ['TARDIS was run only extract Jump Distance (JD) information. The following TARDIS settings were used: Δt bins of '...
        num2str(min(settings.dt_arr)) '-' num2str(max(settings.dt_arr))...
        '; maximum jump distance of ' num2str(settings.maxdist) ' units; '...
        'background frames starting at frame-delay of ' num2str(settings.frame_dist_BG_start) ', using in total ' num2str(size(settings.frame_dist_BG,2)+1) ' frames; '...
        num2str(settings.bgbinningnr) ' BG bins starting at ' num2str(settings.startpointBG) ' units.'];
end

disp(out1)