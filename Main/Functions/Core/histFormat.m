%% General histogram formatting in TARDIS
%---------------------------------------------------------
% Required inputs
% h:            plotted histogram data
%
% Obtained outputs:
% nothing
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function histFormat(h)
    h.EdgeColor = 'k';
    h.EdgeAlpha = 0.2;
    h.FaceColor = [0.5 0.5 0.5];
    h.FaceAlpha = 0.6;
end