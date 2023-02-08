function [J,XB]=index(data,U,center,n,K)
dist = distfcm(center, data);
expo=2;
U=U.^expo;
%% calculate objective function J
J = sum(sum((dist.^2).*U));
%% calculate objective function XB
%%% Calculate denonimator of F
% calculate D1
distc1 = distfcm(center, center);
for i=1:K
distc1(i,i)=999999;
end
D1=min(min(distc1));
D1=D1.^2;
denominator=n*D1;
XB=J/denominator; 
end
function out = distfcm(center, data)
out = zeros(size(center, 1), size(data, 1));
for k = 1:size(center, 1)
    out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end
end



