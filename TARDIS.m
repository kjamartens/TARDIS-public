% -------------------------------------------------------------------------
% Main wrapper for TARDIS (temporal analysis of relative distances in sptPALM)
% Read about TARDIS here:
% https://github.com/kjamartens/TARDIS
% Koen J.A. Martens et al., 2022.
% -------------------------------------------------------------------------
function [parameters, parametersCI, paramEsts, HOfitCI, tottime, time, bgarr, truthoffsetpartial, Visualisation_HO_outputCell,Visualisation_FF_outputCell,anaDDAvisInfoHO,anaDDAvisInfoFF,JDonlydata,SwiftParameters] =...
    TARDIS(poslist,settingsTARDIS)

% Load required folders of functions
addpath(genpath(fileparts(mfilename('fullpath'))));

% Simply run URDA_HO_FF_function.
% Keeping original name for legacy reasons
[parameters, parametersCI, paramEsts, HOfitCI, tottime, time, bgarr, truthoffsetpartial, Visualisation_HO_outputCell,Visualisation_FF_outputCell,anaDDAvisInfoHO,anaDDAvisInfoFF,JDonlydata,SwiftParameters] = ...
    URDA_HO_FF_function(poslist,settingsTARDIS);

end
