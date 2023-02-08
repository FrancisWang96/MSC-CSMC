function [NewS]= SelectMC(X,K,S,c1,c2)
[pY,u,result]=SFC_core_iY(X,K,[],c1,c2);
for t = 1:5
    PC_0(t) = result.validity.PC;
end
num = size(S,1);
flag = zeros(num,1);
for i =1:num
    Si = S(i,:);
    Si = S;
 for t = 1:5
        [pY,u,result]=SFC_core_iY(X,K,Si,c1,c2);
        [indx,~]=valid_clusterIndex(X,pY);
        PC(t) = result.validity.PC;
        CE(t) = result.validity.CE;
 end
        if mean(PC) > mean(PC_0)
            flag(i)=1;
        else flag(i) =0;
        end
       
end
 NewS = S(flag==1,:);
end