%% Replace last line in TARDIS UI
% Replace the last line in TARDIS - usefull for showing fitting progress
% without spamming the UI
%---------------------------------------------------------
% Required inputs
% stringval:    Value that should be showed
% callfromUI:   UI infromation or empty array if displayed generally
%
% Outputs
% None
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function replaceLastTextLine(stringval,callfromUI)
%Check for data in text area - limit to like 200 lines
limitUITextAreaSize(callfromUI,200);
if isempty(callfromUI)
else
    callfromUI.TextArea.Value = callfromUI.TextArea.Value(2:end);
    callfromUI.TextArea.Value = [stringval; callfromUI.TextArea.Value ];
    drawnow
end
end