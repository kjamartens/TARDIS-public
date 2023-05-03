%% Get sizes and input data from JD array
% Get information of DT size, and combine the information in the JD array
%---------------------------------------------------------
% Required inputs
% JDarrSignalCell:      Cell with info about the JDs for every dt
%
% Output
% size_dt:              Size of the dt bins in nr of JDs
% inputarr:             JDs put below one another. Used throughout TARDIS
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [size_dt,inputarr] = getSizesInputFromJDCell(JDarrSignalCell)
size_dt = zeros(1,size(JDarrSignalCell,2));
for dt = 1:size(JDarrSignalCell,2) %Loop over all dt bins wanted
    size_dt(dt) = size(JDarrSignalCell{dt},1);
end
totsize = sum(size_dt);
inputarr = zeros(totsize,1);
for dt = 1:size(JDarrSignalCell,2) %Loop over all dt bins wanted
    if dt == 1
        startpos = 1;
    else
        startpos = sum(size_dt(1:(dt-1)))+1;
    end
    endpos = startpos+size_dt(dt)-1;
    inputarr(startpos:endpos) = JDarrSignalCell{dt}; %Put in a large array combining all dt bins for later use
end