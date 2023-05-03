%% Interpolation of the BG, intra-track linkages - LIN variant
% Function to create an interpolated BG curve. The idea is to get a
% bgbinningnr-by-2 list of x,y positions of the BG curve. Therefore, the
% BGlist is simply histogrammed, and the histogram positions are further
% interpolated makima-style.
%---------------------------------------------------------
% Required inputs
% bgmaxdist:            Maximum distance to lookup BG (m)
% bgbinningnr:          Number of bins used for this
% bginterpolbinsize:    Distnace between bins. Normally
%                       bgmaxdist/bgbinningnr, but sometimes different
% BGlist:               Raw BG linkage data
%
% Output
% output                Interpolated background curve
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [output] = interpolate_BGCurve(bgmaxdist,bgbinningnr,bginterpolbinsize,BGlist)
%Get histogram from BGlist
[h,edges] = histcounts(BGlist,bgbinningnr,'Normalization','probability');
%Get the middle of the bins
midbins = zeros(1,size(edges,2)-1);
for i = 1:size(midbins,2)
    midbins(i) = (edges(i)+edges(i+1))/2;
end
%Interpolate histogram from BGlist to be used later
try
    output(:,2) = [0 interp1(midbins,h,[bginterpolbinsize:bginterpolbinsize:bgmaxdist],'makima')];
    output(:,1) = [0:bginterpolbinsize:bgmaxdist];
catch
    keyboard
end
end