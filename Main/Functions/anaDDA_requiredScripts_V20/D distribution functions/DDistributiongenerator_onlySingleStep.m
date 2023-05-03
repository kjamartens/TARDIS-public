function [framescombined] = DDistributiongenerator_onlySingleStep(koff,kon,Dfree,rangeD,locerrorpdfcorrected,maxindex,fx,fy,maxDindtracking,input,indexframetimerange)
%Probability to start in State 1
 probS1 = kon/(koff+kon);
 %Probability to start in State 2
 probS2 = koff/(kon+koff);
 
frametime = input.frametime;
locerror = input.sigmaerror.^2/frametime;
 
probstationary = exp(-koff*frametime);
probmobile = exp(-kon*frametime);
koff = koff*frametime;
kon = kon*frametime;


%% Integration of 
interpolnumber = 20; %STANDARD THIS IS 10!
x = rangeD(1:interpolnumber:end);
%x = rangeD(1:interpolnumber:maxindex);
x = x';
if input.confinement == 1
    polyfunc = @(x,test) test(9).*x+test(8).*x.^2+test(7).*x.^3+test(6).*x.^4+test(5).*x.^5+test(4).*x.^6+test(3).*x.^7+test(2).*x.^8+test(1).*x.^9+test(10);
    func = @(t1) Integralconfined(t1,Dfree,koff,kon,fx(:,indexframetimerange),fy(:,indexframetimerange),x,polyfunc);
    %func = @(t1) Integralconfined(t1,Dfree,koff,kon,fx,fy,x,polyfunc);
    pdfarray = integral(func,0,1,'ArrayValued',true);
else
    func = @(t1) Integral(t1,Dfree,koff,kon,x,locerror);
    pdfarray = integral(func,0,1,'ArrayValued',true);  
end
% if input.confinement == 1
%     polyfunc = @(x,test) test(5).*x+test(4).*x.^2+test(3).*x.^3+test(2).*x.^4+test(1).*x.^5+test(6);
%     func = @(t1) Integralconfinedtrial(t1,Dfree,koff,kon,fx(:,indexframetimerange),fy(:,indexframetimerange),x,polyfunc);
%     pdfarray = IntegralGaussLobatto(func,0,1,150);
% else
%     %func = @(t1) Integraltrial(t1,Dfree,koff,kon,x,locerror);
%     func = @(t1) IntegraltrialGauss(t1,Dfree,koff,kon,x,locerror);
%     pdfarray = IntegralGaussLobatto(func,0,1,150); 
% end
pdfarray = pdfarray.*(rangeD(2)-rangeD(1));
pdfarray = interp1(pdfarray,1:(1/interpolnumber):length(pdfarray),'spline');

pdfarray(:,2) = pdfarray(:,2)+probstationary.*1/(locerror)*exp(-rangeD(1:length(pdfarray))'./(locerror)).*(rangeD(2)-rangeD(1));
%pdfarray(:,2) = pdfarray(:,2) + probstationary.*normpdf(rangeD(1:length(pdfarray)),0,2*locerror).*(rangeD(2)-rangeD(1));

if input.confinement == 1
pdfarray(:,3) = pdfarray(:,3) + probmobile.*(1./sqrt(polyfunc(Dfree,fx(:,indexframetimerange)))).*(1./sqrt((polyfunc(Dfree,fy(:,indexframetimerange))))).*exp((-(1./polyfunc(Dfree,fx(:,indexframetimerange)))-(1./polyfunc(Dfree,fy(:,indexframetimerange)))).*rangeD(1:length(pdfarray))'/2).*besseli0_fast(((1./polyfunc(Dfree,fx(:,indexframetimerange)))-(1./polyfunc(Dfree,fy(:,indexframetimerange)))).*rangeD(1:length(pdfarray))'/2).*(rangeD(2)-rangeD(1));
%pdfarray(:,3) = pdfarray(:,3) + probmobile.*(1./sqrt(fx(Dfree))).*(1./sqrt(fy(Dfree))).*exp((-(1./fx(Dfree))-(1./fy(Dfree))).*rangeD(1:length(pdfarray))'/2).*besseli0_fast(((1./fx(Dfree)))-(1./fy(Dfree)).*rangeD(1:length(pdfarray))'/2).*(rangeD(2)-rangeD(1));
else
pdfarray(:,3) = pdfarray(:,3) + probmobile.*1/(Dfree+locerror)*exp(-rangeD(1:length(pdfarray))'./(Dfree+locerror)).*(rangeD(2)-rangeD(1));
%pdfarray(:,3) = pdfarray(:,3) + probmobile.*normpdf(rangeD(1:length(pdfarray)),0,2*Dfree+2*locerror).*(rangeD(2)-rangeD(1));

end

% keyboard
% 
% pdfarrayS2S1 = integral(funcodd,0,1,'ArrayValued',true);
% pdfarrayS2S1 = pdfarrayS2S1.*(rangeDStracy(2)-rangeDStracy(1));
% pdfarrayS2S1 = interp1(pdfarrayS2S1,1:(1/interpolnumber):numel(pdfarrayS2S1),'spline')';
% 
% pdfarrayS1S1 = integral(funceven,0,1,'ArrayValued',true);
% pdfarrayS1S1 = pdfarrayS1S1 + probstationary.*1/(locerror)*exp(-x./(locerror));
% pdfarrayS1S1 = pdfarrayS1S1.*(rangeDStracy(2)-rangeDStracy(1));
% pdfarrayS1S1 = interp1(pdfarrayS1S1,1:(1/interpolnumber):numel(pdfarrayS1S1),'spline')';
% 
% pdfarrayS2S2 = integral(funceven2,0,1,'ArrayValued',true);
% 
% if confinement == 1
% pdfarrayS2S2 = pdfarrayS2S2 + probmobile.*(1./sqrt(fx(Dfree))).*(1./sqrt((fy(Dfree)))).*exp((-(1./fx(Dfree))-(1./fy(Dfree))).*x/2).*besseli0_fast(((1./fx(Dfree))-(1./fy(Dfree))).*x/2);
% else
% pdfarrayS2S2 = pdfarrayS2S2 + probmobile.*1/(Dfree+locerror)*exp(-x./(Dfree+locerror));
% end
% pdfarrayS2S2 = pdfarrayS2S2.*(rangeDStracy(2)-rangeDStracy(1));
% pdfarrayS2S2 = interp1(pdfarrayS2S2,1:(1/interpolnumber):numel(pdfarrayS2S2),'spline')';

koff = koff/frametime;
kon = kon/frametime;


if maxDindtracking > 0
    pdfarray(maxDindtracking+1:end,:) = 0;
end

pdfarrayS2S1 = pdfarray(:,1);
pdfarrayS1S1 = pdfarray(:,2);
pdfarrayS2S2 = pdfarray(:,3);    

pdfarrayS1S2 = pdfarrayS2S1*koff/kon;

%% Distribution of D for two frames separated for whether particles starts and finishes within two frames in either state 1 or state 2
Ly=1*(length(pdfarrayS2S2)+length(pdfarrayS2S2))-1;  % 
Ly2=pow2(nextpow2(Ly));    % Find smallest power of 2 that is > Ly

allframes = zeros(maxindex,1);

allframes(1:length(pdfarrayS2S2),1) = probS2*(pdfarrayS2S2+ pdfarrayS2S1) + probS1*(pdfarrayS1S2 +pdfarrayS1S1);

%% Convolution is based on fourier transformation and multiplication
%Fast Fourier transform
% fftpdfarrayS2S2=fft(pdfarrayS2S2, Ly2);		   
% %fftpdfarrayS2S1=fft(pdfarrayS2S1, Ly4);	   % Fast Fourier transform
% fftpdfarrayS1S1=fft(pdfarrayS1S1, Ly2);
% fftpdfarrayS1S2=fft(pdfarrayS1S2, Ly2);
% fftpdfarrayS2S1=kon/koff*fftpdfarrayS1S2;
% 
% %% 2 Frames
% % fftS2S22frames = fftpdfarrayS2S2.*fftpdfarrayS2S2+fftpdfarrayS2S1.*fftpdfarrayS1S2;
% % fftS1S12frames = fftpdfarrayS1S1.*fftpdfarrayS1S1+fftpdfarrayS1S2.*fftpdfarrayS2S1;
% % fftS1S22frames = fftpdfarrayS1S2.*(fftpdfarrayS2S2+fftpdfarrayS1S1);
% % fftS2S12frames = kon/koff*fftS1S22frames;
% StartS12frames = fftpdfarrayS1S1.*fftpdfarrayS1S1+fftpdfarrayS1S2.*(fftpdfarrayS2S2+fftpdfarrayS1S1+fftpdfarrayS2S1);
% StartS22frames = fftpdfarrayS2S2.*fftpdfarrayS2S2+fftpdfarrayS2S1.*(fftpdfarrayS2S2+fftpdfarrayS1S1+fftpdfarrayS1S2);
% 
% fftframescombined = probS2*StartS22frames + probS1* StartS12frames;
% fftframescombined = ifft(fftframescombined, Ly2,'symmetric');
% allframes(:,2) = fftframescombined(1:maxindex);
% 
% %% 3 Frames
% % fftS2S23frames = fftpdfarrayS2S2.*fftS2S22frames+fftpdfarrayS2S1.*fftS1S22frames;
% % fftS1S13frames = fftpdfarrayS1S1.*fftS1S12frames+fftpdfarrayS1S2.*fftS2S12frames;
% % fftS1S23frames = fftpdfarrayS1S2.*fftS2S22frames+fftpdfarrayS1S1.*fftS1S22frames;
% % fftS2S13frames = kon/koff*fftS1S23frames;
% 
% StartS13frames = fftpdfarrayS1S1.*StartS12frames+fftpdfarrayS1S2.*StartS22frames;
% StartS23frames = fftpdfarrayS2S2.*StartS22frames+fftpdfarrayS2S1.*StartS12frames;
%  %fftframescombined = probS2*(fftS2S23frames+ fftS2S13frames) + probS1*(fftS1S23frames +fftS1S13frames);
%  
% fftframescombined = probS1*StartS13frames + probS2*StartS23frames;
% fftframescombined = ifft(fftframescombined, Ly2,'symmetric');
% allframes(:,3) = fftframescombined(1:maxindex);
% 
% %% 4 Frames
% % fftS2S24frames = fftS2S22frames.*fftS2S22frames+  fftS2S12frames.* fftS1S22frames;
% % fftS1S14frames = fftS1S22frames.*fftS2S12frames+  fftS1S12frames.* fftS1S12frames;
% % fftS1S24frames = fftS1S22frames.*(fftS2S22frames+ fftS1S12frames);
% % fftS2S14frames = kon/koff*fftS1S24frames;
%  %StartS14frames =  fftS1S12frames.*fftS1S12frames +fftS1S22frames.*(fftS2S22frames+fftS1S12frames+fftS2S12frames);
%  StartS14frames =  StartS13frames.*fftpdfarrayS1S1+StartS23frames.*fftpdfarrayS1S2;
%  StartS24frames =  StartS23frames.*fftpdfarrayS2S2+StartS13frames.*fftpdfarrayS2S1;
%  %StartS24frames =  fftS2S22frames.*fftS2S22frames +fftS2S12frames.*(fftS2S22frames+fftS1S12frames+fftS1S22frames);
% 
% % combinedstart2 = fftS2S24frames+fftS2S14frames;
% % combinedstart1 = probS1*fftS1S24frames+fftS1S14frames;
% combinedstart2 = probS2*StartS24frames;
% combinedstart1 = probS1*StartS14frames;
% 
% % combinedstart2 = probS2*(fftS2S24frames+fftS2S14frames);
% % combinedstart1 = probS1*(fftS2S14frames+fftS1S14frames);
% 
% fftframescombined = combinedstart2 + combinedstart1;
% fftframescombined = ifft(fftframescombined, Ly2,'symmetric');
% allframes(:,4) = fftframescombined(1:maxindex);
% 
% %% 5 Frames
% % fftS2S25frames = fftS2S22frames.*fftS2S23frames+  fftS2S12frames.* fftS1S23frames;
% % fftS1S15frames = fftS1S22frames.*fftS2S13frames+  fftS1S12frames.* fftS1S13frames;
% % fftS1S25frames = fftS1S22frames.*fftS2S23frames+ fftS1S23frames.*fftS1S12frames;
% % fftS2S15frames = kon/koff*fftS1S25frames;
% % fftframescombined = probS2*(fftS2S25frames+ fftS2S15frames) + probS1*(fftS1S25frames +fftS1S15frames);
% 
% %fftframescombined = (probS2*fftpdfarrayS2S2+probS1*fftpdfarrayS1S2).*combinedstart2+(probS1*fftpdfarrayS1S1+probS2*fftpdfarrayS2S1).*combinedstart1;
% fftframescombined = combinedstart1.*(fftpdfarrayS1S2+fftpdfarrayS1S1)+combinedstart2.*(fftpdfarrayS2S2+fftpdfarrayS2S1);
% fftframescombined = ifft(fftframescombined, Ly2,'symmetric');
% allframes(:,5) = fftframescombined(1:maxindex);
% 
% %% 6 Frames
% % fftS2S26frames = fftS2S23frames.*fftS2S23frames+  fftS2S13frames.* fftS1S23frames;
% % fftS1S16frames = fftS1S23frames.*fftS2S13frames+  fftS1S13frames.* fftS1S13frames;
% % fftS1S26frames = fftS1S23frames.*(fftS2S23frames+ fftS1S13frames);
% % fftS2S16frames = kon/koff*fftS1S26frames;
% % fftframescombined = probS2*(fftS2S26frames+ fftS2S16frames) + probS1*(fftS1S26frames +fftS1S16frames);
% %fftframescombined = (probS2*fftS2S22frames+probS1*fftS1S22frames).*combinedstart2+(probS1*fftS1S12frames+probS2*fftS2S12frames).*combinedstart1;
% fftframescombined = combinedstart1.*StartS12frames+combinedstart2.*StartS22frames;
% fftframescombined = ifft(fftframescombined, Ly2,'symmetric');
% allframes(:,6) = fftframescombined(1:maxindex);
% 
% % fftS2S27frames = fftS2S23frames.*fftS2S24frames+  fftS2S13frames.* fftS1S24frames;
% % fftS1S17frames = fftS1S23frames.*fftS2S14frames+  fftS1S13frames.* fftS1S14frames;
% % fftS1S27frames = fftS1S23frames.*fftS2S24frames+ fftS1S24frames.*fftS1S13frames;
% % fftS2S17frames = kon/koff*fftS1S27frames;
% % fftframescombined = probS2*(fftS2S27frames+ fftS2S17frames) + probS1*(fftS1S27frames +fftS1S17frames);
% %fftframescombined = (probS2*fftS2S23frames+probS1*fftS1S23frames).*combinedstart2+(probS1*fftS1S13frames+probS2*fftS2S13frames).*combinedstart1;
% fftframescombined = combinedstart1.*StartS13frames+combinedstart2.*StartS23frames;
% fftframescombined = ifft(fftframescombined, Ly2,'symmetric');
% allframes(:,7) = fftframescombined(1:maxindex);
% 
% % fftS2S28frames = fftS2S24frames.*fftS2S24frames+  fftS2S14frames.* fftS1S24frames;
% % fftS1S18frames = fftS1S24frames.*fftS2S14frames+  fftS1S14frames.* fftS1S14frames;
% % fftS1S28frames = fftS1S24frames.*(fftS2S24frames+ fftS1S14frames);
% % fftS2S18frames = kon/koff*fftS1S28frames;
% % fftframescombined = probS2*(fftS2S28frames+ fftS2S18frames) + probS1*(fftS1S28frames +fftS1S18frames);
% %fftframescombined = (probS2*fftS2S24frames+probS1*fftS1S24frames).*combinedstart2+(probS1*fftS1S14frames+probS2*fftS2S14frames).*combinedstart1;
% %fftframescombined = combinedstart1.*(fftS1S24frames+fftS1S14frames)+combinedstart2.*(fftS2S24frames+fftS2S14frames);
% fftframescombined = combinedstart1.*StartS14frames+combinedstart2.*StartS24frames;
% fftframescombined = ifft(fftframescombined, Ly2,'symmetric');
% allframes(:,8) = fftframescombined(1:maxindex);

%% Normalize the distributions
%allframes = allframes(1:maxindex,:);%./sum(allframes); 

%% Localization error is slightly differently distributed. This takes this into account by changing distribution 
framerange = 1:1; 
probstationary = probS1*(exp(-koff*framerange*frametime));
allframes = allframes +probstationary.*locerrorpdfcorrected(1:maxindex,1);
%allframes = allframes +probstationary.*locerrorpdfcorrected(1:maxindex,:);            %- probstationary*locerrorpdf;

framescombined = allframes; 
% keyboard   


function [I] = IntegralGaussLobatto(func,lowerbound,upperbound,N)
N = max(3*round((N-1)/3),3) + 1; % Adjust N to the closest valid choice
h = (upperbound - lowerbound)/(N-1);
d = (3/sqrt(5)- 1)*h/2;
t1 = (lowerbound:h:upperbound).'; t1(2:3:N-2) = t1(2:3:N-2) - d; t1(3:3:N-1) = t1(3:3:N-1) + d;
w = ones(1,N); w(4:3:N-3) = 2; w([2:3:N-2,3:3:N-1]) = 5; w = w*h/4;
[A,B,C] = func(t1);
I(:,1) = w*A;
I(:,2) = w*B;
I(:,3) = w*C;
%I = w * func(t1); % Approximately evaluate the integral


function [output] = Integral(t1,Dfree,koff,kon,x,locerror)%,Dfree,koff,kon,fx,fy,x)
output = zeros(size(x,1),3);
dist= 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror));
dist2= 1/((1-t1).*Dfree+locerror)*exp(-x./((1-t1).*Dfree+locerror));
evenconstant = sqrt(koff*kon.*(1-t1)./(t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));

output(:,1) = dist.*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff*t1*(1-t1)));
output(:,2) =dist.*exp(-kon.*t1-koff.*(1-t1)).*evenconstant;
output(:,3) = dist2.*exp(-koff.*t1-kon.*(1-t1)).*evenconstant;


% funcodd = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff*t1*(1-t1)));
% funceven = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-kon.*t1-koff.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));
% funceven2 = @(t1) 1/((1-t1).*Dfree+locerror)*exp(-x./((1-t1).*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-koff.*t1-kon.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));

function [output1, output2, output3] = Integraltrial(t1,Dfree,koff,kon,x,locerror)%,Dfree,koff,kon,fx,fy,x)
%output = zeros(length(t1),length(x),3);
dist= 1./(t1.*Dfree+locerror).*exp(-x./(t1.*Dfree+locerror));
dist2 = flip(dist);
%dist2= 1/((1-t1).*Dfree+locerror)*exp(-x./((1-t1).*Dfree+locerror));
evenconstant = sqrt(koff*kon.*(1-t1)./(t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));
evenconstant(1) = 0;

output1 = dist.*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff.*t1.*(1-t1)));
output2 = dist.*exp(-kon.*t1-koff.*(1-t1)).*evenconstant;
output3 = dist2.*exp(-koff.*t1-kon.*(1-t1)).*evenconstant;
% funcodd = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff*t1*(1-t1)));
% funceven = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-kon.*t1-koff.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));
% funceven2 = @(t1) 1/((1-t1).*Dfree+locerror)*exp(-x./((1-t1).*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-koff.*t1-kon.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));

function [output1, output2, output3] = IntegraltrialGauss(t1,Dfree,koff,kon,x,locerror)%,Dfree,koff,kon,fx,fy,x)
%output = zeros(length(t1),length(x),3);
dist = normpdf(x,0,2*Dfree*t1+2*locerror);
dist2 = normpdf(x,0,2*Dfree*(1-t1)+2*locerror);
%dist= 1./(t1.*Dfree+locerror).*exp(-x./(t1.*Dfree+locerror));
%dist2 = flip(dist);
%dist2= 1/((1-t1).*Dfree+locerror)*exp(-x./((1-t1).*Dfree+locerror));
evenconstant = sqrt(koff*kon.*(1-t1)./(t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));
evenconstant(1) = 0;

output1 = dist.*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff.*t1.*(1-t1)));
output2 = dist.*exp(-kon.*t1-koff.*(1-t1)).*evenconstant;
output3 = dist2.*exp(-koff.*t1-kon.*(1-t1)).*evenconstant;
% funcodd = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff*t1*(1-t1)));
% funceven = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-kon.*t1-koff.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));
% funceven2 = @(t1) 1/((1-t1).*Dfree+locerror)*exp(-x./((1-t1).*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-koff.*t1-kon.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));

function [output1, output2, output3] = Integralconfinedtrial(t1,Dfree,koff,kon,fx,fy,x,polyfunc)%,Dfree,koff,kon,fx,fy,x)
Dx = polyfunc(Dfree.*t1,fx);
Dy = polyfunc(Dfree.*t1,fy);
dist= (1./sqrt(Dx)).*(1./sqrt(Dy)).*exp((-(1./Dx)-(1./Dy)).*x/2).*besseli0_fast(((1./Dx)-(1./Dy)).*x/2);
dist2 = flip(dist,1);
evenconstant = sqrt(koff*kon.*(1-t1)./(t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));
evenconstant(1) = 0;
output1 = dist.*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff.*t1.*(1-t1)));
output2 = dist.*exp(-kon.*t1-koff.*(1-t1)).*evenconstant;
output3 = dist2.*exp(-koff.*t1-kon.*(1-t1)).*evenconstant;
% funcodd = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff*t1*(1-t1)));
% funceven = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-kon.*t1-koff.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));
% funceven2 = @(t1) 1/((1-t1).*Dfree+locerror)*exp(-x./((1-t1).*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-koff.*t1-kon.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));

function [output] = Integralconfined(t1,Dfree,koff,kon,fx,fy,x,polyfunc)
Dx = polyfunc(Dfree.*t1,fx);
Dy = polyfunc(Dfree.*t1,fy);
Dx2 = polyfunc(Dfree.*(1-t1),fx);
Dy2 = polyfunc(Dfree.*(1-t1),fy);
% Dx = fx(Dfree.*t1);
% Dy = fy(Dfree.*t1);
% Dx2 = fx(Dfree.*(1-t1));
% Dy2 = fy(Dfree.*(1-t1));

dist= (1./sqrt(Dx)).*(1./sqrt(Dy)).*exp((-(1./Dx)-(1./Dy)).*x/2).*besseli0_fast(((1./Dx)-(1./Dy)).*x/2);
dist2= (1./sqrt(Dx2)).*(1./sqrt(Dy2)).*exp((-(1./Dx2)-(1./Dy2)).*x/2).*besseli0_fast(((1./Dx2)-(1./Dy2)).*x/2);
evenconstant = sqrt(koff*kon.*(1-t1)./(t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));

output(:,1) = dist.*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff.*t1.*(1-t1)));
output(:,2) =dist.*exp(-kon.*t1-koff.*(1-t1)).*evenconstant;
output(:,3) = dist2.*exp(-koff.*t1-kon.*(1-t1)).*evenconstant;
