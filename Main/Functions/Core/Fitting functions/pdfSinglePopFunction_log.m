%% Generate the likelihood of fitting this pdf function - LOG variant
% The actual fit called by e.g. HistOnly_MLEFit_estimation or
% pdfSinglePopFunction
%---------------------------------------------------------
% Required inputs
% xdata:                All input values of the JD data
% D:                    Diffusion coefficient (m2/s)
% offset:               Offset present in data - important when
%                       BG-subtracted data is not lying at y=0
% loc_unc:              Localization uncertainty in m
% frame_time:           Frame time in s
% strobo_frame_time     Duration that strobo-light is on. If this is
%                       set to -1, the value of the frame will be used instead.
%
% Output
% output                Likelihood of every xdata to belong to these params
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function output = pdfSinglePopFunction_log(xdata,D,offset,loc_unc,strobo_frame_time,frame_time)
% % Add loc uncertainty to D value
% % D = D+loc_unc^2/frame_time;
%Fix strobo frame time in strange conditions
if strobo_frame_time == -1
    strobo_frame_time = frame_time;
elseif strobo_frame_time > frame_time
    disp('Strobo frame time is set higher than frame time! Reverting to frame_time.')
    strobo_frame_time = frame_time;
end

%Following Berglund 2010
% Add loc uncertainty and motion blur coefficient to D value
Dcorr = 4*D*frame_time+2*2*loc_unc^2+4*D*(strobo_frame_time/frame_time)*(1/6)*frame_time;
% Calculate likelihood
output = xdata.*exp(-(xdata.^2)./(Dcorr))+offset;

end