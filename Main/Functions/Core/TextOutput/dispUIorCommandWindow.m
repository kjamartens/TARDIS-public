%% Text visualisation either in UI or in main matlab command window
% Visualisation of strings in UI or main matlab window
%---------------------------------------------------------
% Required inputs
% stringval:    String that should be outputted somewhere
% callfromUI:   UI infromation or empty array if displayed generally

% Outputs
% None
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function dispUIorCommandWindow(stringval,callfromUI)
    %Check for data in text area - limit to like 200 lines
    if isempty(callfromUI)
        disp(stringval)
        global prevMsgLen;
        prevMsgLen = 0;
    else
        limitUITextAreaSize(callfromUI,200);
        callfromUI.TextArea.Value = [stringval; callfromUI.TextArea.Value ];
        drawnow
    end