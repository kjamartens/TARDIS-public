%% Console output replace function
% Output text in the MATLAB console that replaces the existing line
% Used to output fitting values while it's fitting without spamming the
% command window
%---------------------------------------------------------
% Required inputs
% msg:    Message that should be showed now
% Outputs
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function prevMsgLen = consoleOutputReplaceMsg(msg,prevMsgLen)
    global prevMsgLen;
    if isempty(prevMsgLen) 
        prevMsgLen = 0;
    end
    ASCII_BKSP_CHAR = 8;
    if prevMsgLen == 0
        %Check for warnings! Change the previous warning so we can find a
        %new one. This loop only runs at the start of
        %consoleOutputReplaceMsg series.
        [~,warnID] = lastwarn;
        lastwarn('NonsenseID',warnID);
    end
    %Check if a new warning has occured
    if ~strcmp(lastwarn,'NonsenseID')
        %Prevent removal of this warning message
        prevMsgLen = 0;
        [~,warnID] = lastwarn;
        lastwarn('NonsenseID',warnID);
    end
   
    disp([ char(repmat(ASCII_BKSP_CHAR,1,prevMsgLen)) msg]);
    prevMsgLen = numel(msg)+1;
end