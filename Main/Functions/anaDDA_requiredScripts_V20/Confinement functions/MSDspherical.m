function Dx = MSDspherical(r,t,D) 
fun2 = @(x) (x.^2-2).*sin(x)+2*x.*(cos(x));
zerovalues = zeros(100,1);

for i = 2:1:1000
zerovalues(i-1) = fzero(fun2,i);
end

zerovalues = uniquetol(zerovalues,1e-8);
tau = r.^2./D;
zerovalues = zerovalues(zerovalues>0);
start = numel(zerovalues)+1;
zerovalues(start:10000) = zerovalues(end) + pi.*(1:10001-start);
summation = sum(exp(-zerovalues.^2.*t./tau).*1./(zerovalues.^2.*(zerovalues.^2-2)));
MSD = 6/5.*r.^2 - 12.*r.^2.*summation;
Dx = MSD./(6*t);
