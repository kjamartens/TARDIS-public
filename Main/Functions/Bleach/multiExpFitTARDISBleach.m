%Function to get the information from bleach statistics on dt BG bins
%(normalized)
function out = multiExpFitTARDISBleach(lambda,xv,xvnorm)
%     limit = 1000; %limit of 1k is enough, 500 is on the edge -- was 700
%     %initiate bleach values based on this lambda
%     bleachNvalues = 2.^(-(1/(lambda)).*[1:limit]);
%     out = zeros(limit,1);
%     %We loop over the x-positions
%     for n = 1:limit
%         out([1:n-1]) = out([1:n-1])+((n-[1:n-1]).*bleachNvalues([1:n-1]))';
%     end
%     out = out(xv)./max(out(xvnorm));

%     out2 = 2.^(-(1/(lambda*(1+lambda/10)))*(xv+1)); %???????????


    out2 = 2.^(-(1/(lambda*(1)))*(xv+1)); 
    out = ((out2/max(out2)))'; %Normalize




%     keyboard
% Leaving this implementation here because it's understandable, but slower
%     limit = 700; %limit of 1k is enough, 500 is on the edge
%     out = zeros(limit,1);
%     %We loop over the x-positions
%     for n = 1:limit
%         for x = 1:n-1
%             out(x) = out(x)+(n-x)*2^(-(1/lambda)*x);
%         end
%     end
%     out = out(xv)./max(out(xvnorm));
end