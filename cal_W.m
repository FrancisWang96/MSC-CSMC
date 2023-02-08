function [M,C,W] = cal_W(P,Data)
if size(P,1) == 0
    M=[];
    C=[];
    W = zeros(size(Data,1),size(Data,1));
else
    M = P(P(:,3)>0,:);
    C = P(P(:,3)<0,:);
    W = zeros(size(Data,1),size(Data,1));
    pairnum_M= size(M,1);
    pairnum_C= size(C,1);
    Dist = distfcm(Data,Data);
    for i = 1:pairnum_M
        c_ind = int32(M(i,1));
        r_ind = int32(M(i,2));
%         dis = Dist(r_ind,c_ind);
%         w = (2-2*exp(-dis))/(2+2*exp(-dis));
        w=1;
        M(i,3) = w;
        W(c_ind,r_ind) = w;
        W(r_ind,c_ind) = w;
    end
    for i = 1:pairnum_C
        c_ind = C(i,1);
        r_ind = C(i,2);
%         dis = Dist(r_ind,c_ind);
%         w = -4/(2+2*exp(dis));
        w=-1;
        C(i,3) = w;
        W(c_ind,r_ind) = w;
        W(r_ind,c_ind) = w;
    end
end
end
