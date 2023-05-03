%% Change start/upper-lower bound parameters if there's bleach in TARDIS
% Based on orig LB UB
%---------------------------------------------------------
% Required inputs
% bleachStartLbUb:      [start lb ub] of bleach half-time in s
% start:                start params of whatever fitting is being used
% lb:                   lb params of whatever fitting is being used
% ub:                   ub params of whatever fitting is being used
% fittype:              '1pop', '2pop', 'aDDA'
%
% Output
% startparam:           starting parameters
% lowerbound:           lb parameters
% upperbound:           ub parameters
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [startparam, lowerbound, upperbound] = ...
    changeFitParamsBleach(bleachStartLbUb,start,lb,ub,fittype,varargin)
    if strcmp(fittype,'1pop')
        %In this case, start/lb/ub contains [D - BG - BG .... BG]
        %dt = 1 is last entry in start, lb, ub
        startparam = [start(1), start(end), bleachStartLbUb(1)];
        lowerbound = [lb(1), lb(end), bleachStartLbUb(2)];
        upperbound = [ub(1), ub(end), bleachStartLbUb(3)];
    elseif strcmp(fittype,'2pop')
        %In this case, start/lb/ub contains [D1 - D2 - ratio - BG - BG .... BG]
        %dt = 1 is last entry in start, lb, ub
        startparam = [start(1:3), start(end), bleachStartLbUb(1), varargin{1}(1)];
        lowerbound = [lb(1:3), lb(end), bleachStartLbUb(2), varargin{1}(2)];
        upperbound = [ub(1:3), ub(end), bleachStartLbUb(3), varargin{1}(3)];
    elseif strcmp(fittype,'aDDA')
        %In this case, start/lb/ub contains [k1 - k2 - Dfree - BG - BG .... BG]
        %dt = 1 is last entry in start, lb, ub
        startparam = [start(1:3); start(end); bleachStartLbUb(1)];
        lowerbound = [lb(1:3); lb(end); bleachStartLbUb(2)];
        upperbound = [ub(1:3); ub(end); bleachStartLbUb(3)];
    end