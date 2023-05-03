%% Wrapper for aDDA fitting
% Very much taken from Vink et al 2020/2021
%---------------------------------------------------------
% Required inputs

% Obtained outputs:
% output            pdf of aDDA with these input specifications
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [output] = fittinganaDDA(fitspecies, koff, kon, Dfree, rangeD, locerrorpdfcorrected, maxindex, fx, fy, maxDindtracking,input,dt,c,fixedspecies,Dfixed)
combinedpdf = zeros(maxindex,1);
%Loop over the different species
for i = fitspecies
    %Get correct, sigma-error-corrected D1/locerror values
    D1 = input.sigmaerror/dt; 
    locerror = input.sigmaerror*2/dt;
    try
        %Run the aDDA v2.0
        [pdf_allsteps_20]= DDistributiongenerator_v20M(koff(i),kon(i),Dfree(i),D1,rangeD,locerror,fx,fy,maxDindtracking,input,dt);
        pdf = pdf_allsteps_20(:,1);
    catch
        keyboard
    end
    %Normalize the PDF and add it to the other fitspecies
    pdf = pdf./sum(pdf);
    combinedpdf = combinedpdf + c(i)*pdf;
end
%Add fixed species if they exist
if ~isempty(fixedspecies)
    for i = fixedspecies
        combinedpdf = combinedpdf + c(i)*Dfixed(i,1).dist';
    end
end
%Send to output
output = combinedpdf;
end