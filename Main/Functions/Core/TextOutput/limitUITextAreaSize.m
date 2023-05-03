%% Text visualisation either in UI or in main matlab command window
% Function to remove last X lines if there is too much lines there
% I see a pretty tough speed decrease at high nr of lines - observed only at
% like 6000 lines, but I suggest limiting it to 500 or so.
%---------------------------------------------------------
% Required inputs
% callfromUI:   UI infromation or empty array if displayed generally
% nrlines:      Nr of lines to be limited to
% Outputs
% None
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function limitUITextAreaSize(callfromUI,nrlines)
    if max(size(callfromUI.TextArea.Value)) > nrlines
        callfromUI.TextArea.Value(nrlines+1:end) = [];
    end
    drawnow