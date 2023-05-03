function y=Phi2 (bvector,c,xvector,N)
%% Calculation of lauricella series 
% Taken from https://www.semanticscholar.org/paper/A-MATLAB-program-for-the-computation-of-the-%CE%A62-Martos-Naya-Romero-Jerez/2e85394da1207469b9caa8078f6654e2934bbc7b
% and adapted to speed up
 A=15;
% Inverse Laplace Transform
K=exp(A/2);
alphainv =  [0.5,ones(1,N-1)] ;
n =0:N-1;
y1=(-1).^n.*alphainv.*real(LaplacePhi2((A+2*pi*1i*n)/2,c,xvector,bvector));
y=K.*sum(y1,2);
end

function [P]=LaplacePhi2(s,c,x,b)
%P =ones(length(x),length(s));
P = gamma(c)*(s.^(-c));
for k=1:length(b)
P = P./(1-x(:,k)./s);
end
%P = gamma(c)*(s.^(-c)).*P;
end
