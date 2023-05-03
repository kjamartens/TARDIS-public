function [pvalue,KSSTAT]=kstestanaDDA(framenr,parameters,D, input,rangeD,frametimerange)
framescombined = zeros(length(rangeD),framenr);
for i = 1:numel(input.frametimerange)
input.frametime = input.frametimerange(i);
fractionframetimerange = length(D(:,D(3,:)==input.frametime))./length(D);
maxindex = numel(rangeD);
[fx,fy] = Generateconfinedfunction(0:0.05:5,rangeD,input);
fx = fx';
fy = fy';
if input.trackingwindow < 100
      maxD = (input.trackingwindow*input.pixelsize)^2/(4*input.frametime);
      maxDindtracking = round(maxD./(rangeD(2)-rangeD(1)));
else
    maxDindtracking = 0;
end
for ii = 1:size(parameters)
koff = parameters(ii,2);
kon = parameters(ii,3);
Dfree = parameters(ii,4);
c = parameters(ii,1);
framescombinedtemp = DDistributiongenerator(koff,kon,Dfree,rangeD,input.dist(i).locerrorpdfcorrected,maxindex,fx,fy,maxDindtracking,input,1);
framescombinedtemp = framescombinedtemp./sum(framescombinedtemp);
framescombinedtemp = fractionframetimerange*c*framescombinedtemp(:,framenr);
framescombined = framescombined + framescombinedtemp;
end
end
D = D(1,D(2,:)==framenr);
%framescombined = [framescombined(1)/2; framescombined];
func = @(x) interp1(framescombined,x,'spline');
D = sort(D);
cumtrapzframescombined = cumtrapz(framescombined(:,framenr));
cumtrapzframescombined = cumtrapzframescombined + 1-cumtrapzframescombined(end);
Dconverted = (D-rangeD(1))*framenr./(rangeD(3)-rangeD(2))+1;

number =  interp1(cumtrapzframescombined,Dconverted,'spline');
% for i = 1:numel(Dconverted)
% number(i) = integral(func,0.5,Dconverted(i));
% end
number(Dconverted<0.5) =0;
number = sort(number);
try
[ans,pvalue,KSSTAT]=kstest(D,'CDF',[D' number']);
catch
    keyboard
end
% figure
% hold on
% plot(D,number)
% plot(D,1/numel(number):1/numel(number):1)
% figure
% plot(diff([number; 1/numel(number):1/numel(number):1]))

end