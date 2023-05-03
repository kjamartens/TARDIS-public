%% Initialisation of MLE options in TARDIS
%---------------------------------------------------------
% Required inputs
% 
%
% Output
% mleoptions:       Structured MLE options
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function mleoptions = TARDIS_setMLEoptions(invar)
if strcmp(invar,'HO') %Estimation fit
    mleoptions = statset('mlecustom');
    mleoptions.UseParallel = true;
    mleoptions.MaxIter = 1000;
    mleoptions.MaxFunEvals = mleoptions.MaxIter; %Maximum function evaluations
    mleoptions.TolFun = 1e-6; %Termination tolerance on the function value
    mleoptions.TolX = 1e-6; %Termination tolerance for the parameters
    % mleoptions.OptimFun = 'fmincon';
else %Full fit
    mleoptions = statset('mlecustom');
    mleoptions.UseParallel = true;
    mleoptions.MaxIter = 3000;
    mleoptions.MaxFunEvals = mleoptions.MaxIter; %Maximum function evaluations
    mleoptions.TolFun = 1e-12; %Termination tolerance on the function value
    mleoptions.TolX = 1e-9; %Termination tolerance for the parameters
    mleoptions.TolBnd = mleoptions.TolFun*1e-3;
    mleoptions.TolTypeFun = 'abs';
    mleoptions.TolTypeX = 'abs';
end