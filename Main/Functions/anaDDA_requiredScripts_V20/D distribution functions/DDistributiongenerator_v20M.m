function [allframes] = DDistributiongenerator_v20M(koff,kon,Dfree,D1,rangeD,locerror,fx,fy,maxDindtracking,input,indexframetimerange)
% keyboard
%%

%%
frametime = input.frametime;
koff = koff*frametime;
kon = kon*frametime;

%% Integration of 
interval = (log(rangeD(end))-log(rangeD(1)))/input.integrationinterval; 
x = exp(log(rangeD(1))-2*interval:interval:log(rangeD(end))+2*interval);
x = x';
if input.confinement == 1
    polyfunc = @(x,test) test(9).*x+test(8).*x.^2+test(7).*x.^3+test(6).*x.^4+test(5).*x.^5+test(4).*x.^6+test(3).*x.^7+test(2).*x.^8+test(1).*x.^9+test(10);
    %setting all fx to 1
    func = @(t1) Integralconfined(t1,Dfree,locerror,koff,kon,fx(:,indexframetimerange),fy(:,indexframetimerange),x,polyfunc);
%     func = @(t1) Integralconfined(t1,Dfree,locerror,koff,kon,fx(:,1),fy(:,1),x,polyfunc);
    pdfarray = integral(func,0,1,'ArrayValued',true);
else
    func = @(t1) Integral(t1,Dfree,D1,koff,kon,x,locerror);
    pdfarray = integral(func,0,1,'ArrayValued',true); 
    polyfunc = NaN;
end

pdfarray = pdfarray.*(rangeD(2)-rangeD(1));
% 
% rangeDprecise = rangeD(1):(rangeD(2)-rangeD(1)):rangeD(end)/100;
% pdfarray2 = interp1(x,pdfarray,rangeDprecise,'spline');
% allframesprecise = convolutedistributions(koff,kon,Dfree,D1,locerror,pdfarray2,rangeDprecise,fx,fy,maxDindtracking,input,indexframetimerange);
% rangeDfast = rangeD(1)*10:(rangeD(2)-rangeD(1))*10:rangeD(end)/10+2*((rangeD(2)-rangeD(1))*10);
% pdfarray2 = interp1(x,pdfarray,rangeDfast,'spline');
% allframesprecise = convolutedistributions(koff,kon,Dfree,D1,locerror,pdfarray2,rangeDfast,fx,fy,maxDindtracking,input,indexframetimerange);
% rangeDfast = rangeD(1)*100:(rangeD(2)-rangeD(1))*100:rangeD(end)+2*((rangeD(2)-rangeD(1))*100);
% pdfarray2 = interp1(x,pdfarray,rangeDfast,'spline')*10;
% allframes = convolutedistributions(koff,kon,Dfree,D1,locerror,pdfarray2,rangeDfast,fx,fy,maxDindtracking,input,indexframetimerange);
% rangeDinterp = (rangeDprecise(end)+(rangeD(2)-rangeD(1))):(rangeD(2)-rangeD(1)):rangeD(end);
% allframes = [allframesprecise; interp1(rangeDfast,allframes,rangeDinterp)/10];

pdfarray = interp1(x,pdfarray,rangeD,'spline');
allframes = convolutedistributions(koff,kon,Dfree,D1,locerror,pdfarray,rangeD,fx,fy,maxDindtracking,input,indexframetimerange,polyfunc);

function allframes = convolutedistributions(koff,kon,Dfree,D1,locerror,pdfarray,rangeD,fx,fy,maxDindtracking,input,indexframetimerange,polyfunc)
frametime = input.frametime;
probstationary = exp(-koff);
probmobile = exp(-kon);

%Probability to start in State 1
probS1 = kon/(koff+kon);
%Probability to start in State 2
probS2 = koff/(kon+koff);
pdfarray(:,2) = pdfarray(:,2)+probstationary.*1/(D1+locerror)*exp(-rangeD'./(D1+locerror)).*(rangeD(2)-rangeD(1));
%pdfarray(:,2) = pdfarray(:,2) + probstationary.*normpdf(rangeD(1:length(pdfarray)),0,2*locerror).*(rangeD(2)-rangeD(1));

if input.confinement == 1
pdfarray(:,3) = pdfarray(:,3) + probmobile.*(1./sqrt(polyfunc(Dfree,fx(:,indexframetimerange)))).*(1./sqrt((polyfunc(Dfree,fy(:,indexframetimerange))))).*exp((-(1./polyfunc(Dfree,fx(:,indexframetimerange)))-(1./polyfunc(Dfree,fy(:,indexframetimerange)))).*rangeD(1:length(pdfarray))'/2).*besseli0_fast(((1./polyfunc(Dfree,fx(:,indexframetimerange)))-(1./polyfunc(Dfree,fy(:,indexframetimerange)))).*rangeD(1:length(pdfarray))'/2).*(rangeD(2)-rangeD(1));
%pdfarray(:,3) = pdfarray(:,3) + probmobile.*(1./sqrt(fx(Dfree))).*(1./sqrt(fy(Dfree))).*exp((-(1./fx(Dfree))-(1./fy(Dfree))).*rangeD(1:length(pdfarray))'/2).*besseli0_fast(((1./fx(Dfree)))-(1./fy(Dfree)).*rangeD(1:length(pdfarray))'/2).*(rangeD(2)-rangeD(1));
else
pdfarray(:,3) = pdfarray(:,3) + probmobile.*1/(Dfree+locerror)*exp(-rangeD'./(Dfree+locerror)).*(rangeD(2)-rangeD(1));
%pdfarray(:,3) = pdfarray(:,3) + probmobile.*normpdf(rangeD(1:length(pdfarray)),0,2*Dfree+2*locerror).*(rangeD(2)-rangeD(1));
end

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
combdist = probS2*(pdfarrayS2S2+ pdfarrayS2S1) + probS1*(pdfarrayS1S2 +pdfarrayS1S1);

maxindex = length(rangeD);

%Ly2=pow2(nextpow2(Ly-1)); % Find smallest power of 2 that is => Ly
allframes = zeros(maxindex,8);

%allframes(1:length(pdfarrayS2S2),1) = probS2*(pdfarrayS2S2+ pdfarrayS2S1) + probS1*(pdfarrayS1S2 +pdfarrayS1S1);
allframes(:,1) = combdist;

%% Convolution is based on fourier transformation and multiplication
%Fast Fourier transform
% fftpdfarrayS2S2=fft(pdfarrayS2S2, maxindex);		   
% %fftpdfarrayS2S1=fft(pdfarrayS2S1, Ly4);	   % Fast Fourier transform
% fftpdfarrayS1S1=fft(pdfarrayS1S1, maxindex);
% fftpdfarrayS1S2=fft(pdfarrayS1S2, maxindex);
% fftpdfarrayS2S1=kon/koff*fftpdfarrayS1S2;

% %% 2 Frames
% % fftS2S22frames = fftpdfarrayS2S2.*fftpdfarrayS2S2+fftpdfarrayS2S1.*fftpdfarrayS1S2;
% % fftS1S12frames = fftpdfarrayS1S1.*fftpdfarrayS1S1+fftpdfarrayS1S2.*fftpdfarrayS2S1;
% % fftS1S22frames = fftpdfarrayS1S2.*(fftpdfarrayS2S2+fftpdfarrayS1S1);
% % fftS2S12frames = kon/koff*fftS1S22frames;
% StartS12frames = fftpdfarrayS1S1.*fftpdfarrayS1S1+fftpdfarrayS1S2.*(fftpdfarrayS2S2+fftpdfarrayS1S1+fftpdfarrayS2S1);
% StartS22frames = fftpdfarrayS2S2.*fftpdfarrayS2S2+fftpdfarrayS2S1.*(fftpdfarrayS2S2+fftpdfarrayS1S1+fftpdfarrayS1S2);
% 
% fftframescombined = probS2*StartS22frames + probS1* StartS12frames;
% fftframescombined = ifft(fftframescombined, maxindex,'symmetric');
% allframes(:,2) = fftframescombined;
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
% fftframescombined = ifft(fftframescombined, maxindex,'symmetric');
% allframes(:,3) = fftframescombined;
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
% fftframescombined = ifft(fftframescombined, maxindex,'symmetric');
% allframes(:,4) = fftframescombined;
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
% fftframescombined = ifft(fftframescombined, maxindex,'symmetric');
% allframes(:,5) = fftframescombined;
% 
% %% 6 Frames
% % fftS2S26frames = fftS2S23frames.*fftS2S23frames+  fftS2S13frames.* fftS1S23frames;
% % fftS1S16frames = fftS1S23frames.*fftS2S13frames+  fftS1S13frames.* fftS1S13frames;
% % fftS1S26frames = fftS1S23frames.*(fftS2S23frames+ fftS1S13frames);
% % fftS2S16frames = kon/koff*fftS1S26frames;
% % fftframescombined = probS2*(fftS2S26frames+ fftS2S16frames) + probS1*(fftS1S26frames +fftS1S16frames);
% %fftframescombined = (probS2*fftS2S22frames+probS1*fftS1S22frames).*combinedstart2+(probS1*fftS1S12frames+probS2*fftS2S12frames).*combinedstart1;
% fftframescombined = combinedstart1.*StartS12frames+combinedstart2.*StartS22frames;
% fftframescombined = ifft(fftframescombined, maxindex,'symmetric');
% allframes(:,6) = fftframescombined;
% 
% % fftS2S27frames = fftS2S23frames.*fftS2S24frames+  fftS2S13frames.* fftS1S24frames;
% % fftS1S17frames = fftS1S23frames.*fftS2S14frames+  fftS1S13frames.* fftS1S14frames;
% % fftS1S27frames = fftS1S23frames.*fftS2S24frames+ fftS1S24frames.*fftS1S13frames;
% % fftS2S17frames = kon/koff*fftS1S27frames;
% % fftframescombined = probS2*(fftS2S27frames+ fftS2S17frames) + probS1*(fftS1S27frames +fftS1S17frames);
% %fftframescombined = (probS2*fftS2S23frames+probS1*fftS1S23frames).*combinedstart2+(probS1*fftS1S13frames+probS2*fftS2S13frames).*combinedstart1;
% fftframescombined = combinedstart1.*StartS13frames+combinedstart2.*StartS23frames;
% fftframescombined = ifft(fftframescombined, maxindex,'symmetric');
% allframes(:,7) = fftframescombined;
% 
% % fftS2S28frames = fftS2S24frames.*fftS2S24frames+  fftS2S14frames.* fftS1S24frames;
% % fftS1S18frames = fftS1S24frames.*fftS2S14frames+  fftS1S14frames.* fftS1S14frames;
% % fftS1S28frames = fftS1S24frames.*(fftS2S24frames+ fftS1S14frames);
% % fftS2S18frames = kon/koff*fftS1S28frames;
% % fftframescombined = probS2*(fftS2S28frames+ fftS2S18frames) + probS1*(fftS1S28frames +fftS1S18frames);
% %fftframescombined = (probS2*fftS2S24frames+probS1*fftS1S24frames).*combinedstart2+(probS1*fftS1S14frames+probS2*fftS2S14frames).*combinedstart1;
% %fftframescombined = combinedstart1.*(fftS1S24frames+fftS1S14frames)+combinedstart2.*(fftS2S24frames+fftS2S14frames);
% fftframescombined = combinedstart1.*StartS14frames+combinedstart2.*StartS24frames;
% fftframescombined = ifft(fftframescombined, maxindex,'symmetric');
% allframes(:,8) = fftframescombined;

%% Normalize the distributions
framerange = 2:8; 
%allframes = allframes(1:maxindex,:);%./sum(allframes); 
rangex = input.rangex;
locdist = input.locdist;
conversion = locerror*(rangex(2,2)-rangex(1,2))/(2*(rangeD(2)-rangeD(1)));
maxindexloc = find(rangeD/locerror>max(rangex(:)),1)-1;
for i = 8:-1:2
    %corrected(1:maxindexloc,i)=interp1(rangex(:,i),locdist(:,i),rangeD(1:maxindexloc)/locerror')/(conversion*i);
    corrected(1:maxindexloc,i)=locdist{i}(rangeD(1:maxindexloc)/locerror')/(conversion*i);
end
corrected(isnan(corrected))=0;
locerrorpdf = gammapdf(rangeD(1:length(corrected))',2:8,(D1+locerror)).*(rangeD(2)-rangeD(1));
locerrorpdf(corrected(:,2:8)==0)=0;
locerrorpdfcorrected = corrected(:,2:8)-locerrorpdf;
%% Localization error is slightly differently distributed. This takes this into account by changing distribution 
probstationary = probS1*(exp(-koff*framerange*frametime));
if D1 == 0
    allframes(1:length(corrected),2:8) = allframes(1:length(corrected),2:8)+probstationary.*locerrorpdfcorrected;
end
allframes(maxindex+1:numel(rangeD),:)=1e-25;

function [output] = Integral(t1,Dfree,D1,koff,kon,x,locerror)%,Dfree,koff,kon,fx,fy,x)
%output = zeros(numel(x),3);
output= 1/(t1.*Dfree+(1-t1).*D1+locerror)*exp(-x./(t1.*Dfree+(1-t1).*D1+locerror)).*exp(-kon.*t1-koff.*(1-t1));
%dist2= 1/((1-t1).*Dfree+t1.*D1+locerror)*exp(-x./((1-t1).*Dfree+t1.*D1+locerror));
evenconstant = sqrt(koff*kon.*(1-t1)./(t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));

output = [output output.*evenconstant];
output(:,1) = output(:,1).*kon.*besseli(0,2*sqrt(kon*koff*t1*(1-t1)));

%output(:,3) = dist2.*exp(-koff.*t1-kon.*(1-t1)).*evenconstant;
output = [output 1/((1-t1).*Dfree+t1.*D1+locerror)*exp(-x./((1-t1).*Dfree+t1.*D1+locerror)).*exp(-koff.*t1-kon.*(1-t1)).*evenconstant];
% funcodd = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*kon.*exp(-kon.*t1-koff.*(1-t1)).*besseli(0,2*sqrt(kon*koff*t1*(1-t1)));
% funceven = @(t1) 1/(t1.*Dfree+locerror)*exp(-x./(t1.*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-kon.*t1-koff.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));
% funceven2 = @(t1) 1/((1-t1).*Dfree+locerror)*exp(-x./((1-t1).*Dfree+locerror)).*sqrt(koff*kon.*(1-t1)./(t1)).*exp(-koff.*t1-kon.*(1-t1)).*besseli(1,2*sqrt(kon*koff.*t1.*(1-t1)));

function [output] = Integralconfined(t1,Dfree,locerror,koff,kon,fx,fy,x,polyfunc)
Dx = polyfunc(Dfree.*t1,fx)+locerror;
Dy = polyfunc(Dfree.*t1,fy)+locerror;
Dx2 = polyfunc(Dfree.*(1-t1),fx)+locerror;
Dy2 = polyfunc(Dfree.*(1-t1),fy)+locerror;
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

function dist = gammapdf(x,a,b)
dist = 1./(b.^a.*gamma(a)).*x.^(a-1).*exp(-x./b);