%% Rounding digit numbers
% A function to ensure the same number of characters are used for a value
% representation. Leading values will be spaces, trailing values will be 0s
%---------------------------------------------------------
% Required inputs
% value:                Value to be rounded
% nrdigitsprepoint:     Nr of digits before the period
% nrdigitspastpoint:    Nr of digits after the period

% Outputs
% string:               Properly formatted string of the value
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function string = roundNrDigits(value,nrdigitsprepoint,nrdigitspastpoint)
%Weird floating error fix
value = value+1e-9;
%Find the current number of digits before decimal point
nrdigitsprepointcurr = length(num2str(floor(value)));
%Initiate trailing info
leadinfo = '';
if nrdigitsprepointcurr < nrdigitsprepoint
    for i = 1:(nrdigitsprepoint-nrdigitsprepointcurr)
        leadinfo = [leadinfo ' '];
    end
end
%Get trailing info
inputtrailstring = num2str(mod(round(value,nrdigitspastpoint),1));
inputtrailstring = inputtrailstring(3:end);
if length(inputtrailstring)<nrdigitspastpoint
    for i = 1:(nrdigitspastpoint-length(inputtrailstring))
        inputtrailstring = [inputtrailstring '0'];
    end
end
%Create final string

%Pre-dot value needs to be increased by half of trailing digits
if nrdigitspastpoint > 0
    string = [leadinfo num2str(floor(value+5/10^(1+nrdigitspastpoint))) '.' inputtrailstring];
else
    string = [leadinfo num2str(floor(value+5/10^(1+nrdigitspastpoint)))];
end
end