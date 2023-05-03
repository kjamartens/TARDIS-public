
Weighted Histogram Utility

 Mehmet Suzen
 mehmet dot suzen at physics dot org
 (c) 2013, BSD license

histwc is a utility to construct weighted histogram. Given
number of bins, values and corresponding weights, there are 
two functions produces identical outputs: a vector of histogram 
frequencies and intervals.  Functions generate a vector of 
cumulative weights for data histogram. Equal number of bins 
will be considered using minimum and maximum values of the data. 
Weights will be summed in the given bin. These two vectors can 
be used to plot the weighted histogram. One of the function is 
vectorized, but not necessarily faster.

Example:
%vv    = [1.0 1.0 2.0 5.0 3.0 4.0 1.0 5.0 3.0 1.0];     % values
%ww    = [0.1 0.1 0.05 0.05 0.2 0.05 0.05 0.1 0.1 0.2]; % weights
%nbins = 5;
%[histw, intervals] = histwc(vv, ww, nbins);
% Example Visualise:
% bar(intervals, histw)
