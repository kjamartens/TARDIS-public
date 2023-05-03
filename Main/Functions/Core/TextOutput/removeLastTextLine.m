%% Remove last line in UI
% remove the last line in the TARDIS UI
%---------------------------------------------------------
% Required inputs
% callfromUI:    UI info of TARDIS
% Outputs
% None
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function removeLastTextLine(callfromUI)
if isempty(callfromUI)
else
    callfromUI.TextArea.Value = callfromUI.TextArea.Value(2:end);
    drawnow
end
end