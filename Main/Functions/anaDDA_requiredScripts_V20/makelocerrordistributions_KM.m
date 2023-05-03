function [locerrorpdf,locerrorpdfcorrected] = makelocerrordistributions_KM(rangeD,locerror,input)
% This function produces the distributions of an uncorrected and corrected
% immobile particles that generate an apparent mobility due to localization
% error. The uncorrected pdf is later substracted from the total
% distribution and replaced by the corrected pdf. 
% Inputs:
% rangeD: is the range of x values for which the distribution is calculated
% locerror: is the localization error^2/dt (in um^2/s)
% Outputs:
% locerrorpdf: The uncorrected distribution 
% locerrorpdfcorrected: the corrected distribution
% Created 07-08-2019 by Jochem Vink

%maxcalculatedD = 2;
Dsigmaerror = input.sigmaerror^2/min(input.frametimerange);
maxcalculatedD =  -log(Dsigmaerror*1e-10)*Dsigmaerror;
maxrange = numel(rangeD(1)/max(input.framerange):((rangeD(2)-rangeD(1))/max(input.framerange)):max(rangeD(:)));
locerrorpdf = zeros(maxrange,8);
locerrorpdfcorrected = zeros(maxrange,8);

for i = input.framerange
    binsize = (rangeD(2)-rangeD(1))/i;     % Binsize between values of D to be calculated is adjusted for amount of steps of data 
    adjustedrangeD = rangeD(1)/i:binsize:max(rangeD(:)); 
    adjustedrangeD = adjustedrangeD';
    
    locerrorpdf(1:numel(adjustedrangeD),i) = ((i./locerror).^i).*adjustedrangeD.^(i-1).*(exp(-(i.*adjustedrangeD/locerror))./gamma(i)).*binsize;
    locerrorpdf(end:maxrange,i) = 0;
    precisionfactor = 1;
    pass = 0;
    if i > 1  % For i = 1 steps are uncorrelated so no correction required, for i>1 lauricella series are required for calculation of distribution
        while pass == 0
        range = rangeD(1)/i:binsize:maxcalculatedD; 
        rangefit = range(1:precisionfactor:end);     % To speed up calculation only calculate part of the distribution and calculate others via interpolation
        
        corrected=laur(rangefit,1,locerror/i,i).*binsize;    % Calculate values with lauricella series
        corrected = interp1(corrected,1:1/precisionfactor:numel(corrected));
        
        if i == 2
        if sum(corrected)>1.0005
            warning('With current parameters error is relatively large: try increasing precision parameter in the input file')
%             precision = inputdlg('With current parameters error is relatively large: try increasing precision','Precision increase',1,{num2str(input.precision)})
%             prompt(['precision increased from ' num2str(input.precision) ' to ' precision])
%             input.precision = str2num(
        end    
        if sum(corrected)<0.9995
            maxcalculatedD = maxcalculatedD*2;
        else
            pass = 1;
        end
        else
            pass = 1;
        end
        
        locerrorpdfcorrected(1:numel(corrected),i)=corrected;
        locerrorpdf(numel(corrected)+1:numel(adjustedrangeD),i)=0;
        end
    else
        locerrorpdfcorrected(1:numel(adjustedrangeD),i) =((i./locerror).^i).*adjustedrangeD.^(i-1).*(exp(-(i.*adjustedrangeD/locerror))./gamma(i)).*binsize;
    end
  
end
locerrorpdfcorrected(:,input.framerange) = locerrorpdfcorrected(:,input.framerange)./factorial(input.framerange-1);
locerrorpdfcorrected = locerrorpdfcorrected./sum(locerrorpdfcorrected);
locerrorpdf = locerrorpdf./sum(locerrorpdf);
locerrorpdfcorrected = locerrorpdfcorrected-locerrorpdf;