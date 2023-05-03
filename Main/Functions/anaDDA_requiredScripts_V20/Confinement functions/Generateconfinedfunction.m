function [fx,fy] = Generateconfinedfunction(rangeDPDA,rangeDStracy,input)
% This function calculates the distribution of D given a certain cell
% geometry, localization error and frametime. This distribution is derived from equations 
t = input.frametime;
locerror = input.sigmaerror.^2;

%% Generate MSDs for X confinement
% Size of box in x dimension is calculated
xbox = sqrt(pi()*input.radiusofcell.^2);        

a = xbox.^2/6;
b = 16*xbox.^2/(pi^4);
c = 0;

for n = 1:2:5905
tempc = (1./n.^4)*exp(-0.5*(n*pi./xbox).^2.*(rangeDPDA*2)*t);
c = c+tempc;
end

MSDx = a-b*c+2*locerror;
Dx = MSDx./(2*t);
fun = @(x) 0.5*(besselj(0,x)-besselj(2,x));

for i = 1:100
zerovalues(i) = fzero(fun,i);
end

zerovalues = uniquetol(zerovalues,1e-8);
Volumecylinder = pi().*input.radiusofcell.^2.*(input.lengthcell-2.*input.radiusofcell);
Volumesphere = 4./3*pi().*input.radiusofcell.^3;
Dx = MSDcircular(zerovalues,input.radiusofcell,t,rangeDPDA)./(4.*t);
Dx = Dx.*Volumecylinder./(Volumecylinder+Volumesphere);
fun = @(r) MSDcircular(zerovalues,sqrt(input.radiusofcell.^2-r.^2),t,rangeDPDA).*pi().*(input.radiusofcell.^2-r.^2);
%Dxsphere = integral(fun,-input.radiusofcell,input.radiusofcell,'Arrayvalued',true)./(4.*t);
Dxsphere = MSDspherical(input.radiusofcell,t,rangeDPDA);
Dxsphere = Dxsphere.*Volumesphere./(Volumecylinder+Volumesphere);
Dx = Dxsphere + Dx';
%Dx = Dx';
Dx = Dx + input.sigmaerror.^2./input.frametime;


ybox = input.lengthcell-2*input.radiusofcell+input.radiusofcell*(4/3*pi).^(1/3);

a = ybox.^2/6;
b = 16*ybox.^2/(pi^4);

c = 0;

for n = 1:2:10905
tempc = (1/n.^4)*exp(-0.5*(n*pi./ybox).^2.*(rangeDPDA*2)*t);
c = c+tempc;
end

MSDy = a-b*c+2*locerror;
Dy = MSDy./(2*t);
% Dy2 = rangeDPDA.*Volumecylinder./(Volumecylinder+Volumesphere);
% Dy2 = Dy2 + Dxsphere;
% Dy2 = Dy2 + input.sigmaerror.^2/input.frametime;
% Dy = Dy2;
if abs(input.lengthcell - 2.*input.radiusofcell) <0.1
Dy = Dx;
end

% %fun = @(r) MSDbox((input.lengthcell-2*input.radiusofcell)+2*sqrt(input.radiusofcell.^2-r.^2),rangeDPDA,0.01).*((input.lengthcell-2*input.radiusofcell)+2*sqrt(input.radiusofcell.^2-r.^2)).*2*pi*r;
% fun = @(r) MSDbox(input.lengthcell-2*input.radiusofcell+2*r,rangeDPDA,0.01).*pi().*(input.radiusofcell.^2-r.^2);
% Dy2 = integral(fun,-input.radiusofcell,input.radiusofcell,'Arrayvalued',true);
% %Dy2 = Dy2./(Volumecylinder+Volumesphere);
% Dy2 = Dy2./Volumesphere;
% Dy2 = Dy2 + input.sigmaerror.^2/input.frametime;
% Dy = rangeDPDA.*Volumecylinder./(Volumecylinder+Volumesphere);
% Dy = Dy + Dxsphere'*Volumesphere./(Volumecylinder+Volumesphere);
% Dy = Dy+ input.sigmaerror.^2/input.frametime;

%Dy = Dy2;
x = rangeDStracy';
% pdfarray = (1./sqrt(Dx)).*(1./sqrt(Dy)).*exp((-(1./Dx)-(1./Dy)).*x/2).*besseli0_fast(((1./Dx)-(1./Dy)).*x/2)*(x(2)-x(1));
% pdfarray = single(pdfarray);
% pdfarray = gpuArray(pdfarray);
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
fx = polyfit(rangeDPDA(2:end)',Dx(2:end)',9);
fy = polyfit(rangeDPDA(2:end)',Dy(2:end)',9);
% fx = fit(rangeDPDA(2:end)',Dx(2:end)','poly9');
% fy = fit(rangeDPDA(2:end)',Dy(2:end)','poly9');
% fx = fit(rangeDPDA',Dx','smoothingspline');
% fy = fit(rangeDPDA',Dy','smoothingspline');
% for i = 1:length(MSDy)
% pdfarray(:,i) = fconv(func(x,MSDx(i)),func(x,MSDy(i)));
% %pdfarray(:,i) = fconv(gampdf([rangeDStracy(2)*t 2*t*rangeDStracy(2:length(rangeDStracy))],0.5,2*MSDx(i)),gampdf([t*rangeDStracy(2) 2*t*rangeDStracy(2:length(rangeDStracy))],0.5,2*MSDy(i)));
% end
%fx = [fx.p1 fx.p2 fx.p3 fx.p4 fx.p5 fx.p6 fx.p7 fx.p8 fx.p9 fx.p10];
%fy = [fy.p1 fy.p2 fy.p3 fy.p4 fy.p5 fy.p6 fy.p7 fy.p8 fy.p9 fy.p10];
