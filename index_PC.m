function [J,XB,c]=index_PC(ipop,u,tc,K)
S=ipop.pairS;
X=ipop.data;
c1=ipop.weight(1);
c2=ipop.weight(2);
dist=distfcm(tc,X);
mf = u.^2;
sumc = sum(sum((dist.^2-c2).*mf));
jm1=sumc;
c=0;
for i=1:size(S,1)
    up = S(i,1);
    uq = S(i,2);
    sval = S(i,3);
    if sval>=0
       c =c+sum(0.5*sval*(u(:,up)-u(:,uq)).^2);
    else
       c = c+sum(-1.0*sval*dot(u(:,up),u(:,uq)));
    end
end
jm2=c*c1;
J=jm1+c*c1;  
%% calculate objective function XB
%%% Calculate denonimator of F
% calculate D1
distc1 = distfcm(tc,tc).^2;
for i=1:K
    distc1(i,i)=999999;
end
D1=min(min(distc1));
denominator=size(X,1)*D1;
XB=jm1/denominator; 
end