%% Interpolation of the BG, intra-track linkages - LOG variant
% Function to create an interpolated BG curve. The idea is to get a
% bgbinningnr-by-2 list of x,y positions of the BG curve. Therefore, the
% BGlist is simply histogrammed, and the histogram positions are further
% interpolated makima-style.
%---------------------------------------------------------
% Required inputs
% bgmaxdist:            Maximum distance to lookup BG (m)
% bgmindist:            Minimum distance to lookup BG (m)
% bgbinningnr:          Number of bins used for this
% bgbinningnr_interpolate:   size of output
% BGlist:               Raw BG linkage data
%
% Output
% output                Interpolated background curve
%---------------------------------------------------------
% Koen J.A. Martens, 2022
%---------------------------------------------------------
function [output] = interpolate_BGCurve_log(bgmaxdist,bgmindist,bgbinningnr,bgbinningnr_interpolate,BGlist)
%Get histogram from BGlist
[h,edges] = histcounts(BGlist,logspace(log10(bgmindist),log10(bgmaxdist),bgbinningnr+1),'Normalization','probability');
% [h,edges] = histcounts(BGlist,bgbinningnr,'Normalization','probability');
%Get the middle of the bins
midbins = zeros(1,size(edges,2)-1);
for i = 1:size(midbins,2)
    midbins(i) = sqrt(edges(i).*edges(i+1));
end
%Interpolate histogram from BGlist to be used later
try
    output(:,2) = interp1(midbins,h,logspace(log10(bgmindist),log10(bgmaxdist),bgbinningnr_interpolate+1),'makima');
    output(:,1) = logspace(log10(bgmindist),log10(bgmaxdist),bgbinningnr_interpolate+1);
catch
    keyboard
end
end