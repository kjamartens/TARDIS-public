function out = Olkin_Pratt_Estimator(R2,N,p,K)
%% Determines the estimated O-P statistical value
% Beginning part - standard
p1 = (N-3)/(N-p-1)*(1-R2);

% Hypergeometic part
a = 1;
b = 1;
c = (N-p+1)/2;
z = 1-R2;

%Calculate tk over k = 1..K
for k = 0:K
    if k > 0
        %Rising factorial
        ak = a;
        bk = b;
        ck = c;
        for kk = 1:(k-1)
            ak = ak*(a+kk);
            bk = bk*(b+kk);
            ck = ck*(c+kk);
        end
    else
        ak = 1;
        bk = 1;
        ck = 1;
    end
    %Calculation of t_k
    t(k+1) = ((ak*bk)/ck)*((z^k)/(factorial(k)));

    %Here also calculating the formula by Shieh, but I stick with the top
    %one - differences are incredibly small (order of 1e-8) (but this
    %method is slightly faster, which we don't need)
%     if k == 0
%         t2(k+1) = 1;
%     else
%         t2(k+1) = ((gamma(1+k)*gamma((N-p-1)/2))/gamma((N-p-1)/2+k))*(1-R2)^k;
%     end
end
hygeo = sum(t);
out = 1-p1*hygeo;

% Code for the exact Olkinn-Pratt. Works, but not used because of wrong floating-point
% arithmatic at R2 > 0.5
%Also at R2 < 0.5, differences are in the order of 1e-12

% if z == 0
%     exactOP = 1;
% elseif z == 1
%     exactOP = (c-1)/(c-2);
% elseif mod(c,1) == 0
%     shared = ((z-1)/z);
%     the_sum = 0;
%     for k = 2:(c-1)
%         the_sum = the_sum + (shared^k)/(c-k);
%     end
%     prefactor = (c-1)*z*(z-1)^(-2);
%     exactOP = prefactor*(the_sum-shared^c*log(1-z));
% else
%     cur_res = 1/(1-z)*(1+(sqrt(z)*asin(sqrt(z)))/sqrt(1-z));
%     for i = 1.5:1:c
%         cur_c = i-1;
%         cur_res = (cur_c-cur_c*(1-z)*cur_res)/(z*(cur_c-1));
%     end
%     exactOP = cur_res;
% end
end