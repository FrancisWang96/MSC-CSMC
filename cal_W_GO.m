function [M,C,W] = cal_W_GO(P,Data,GOsim,xita)
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
        sim = GOsim(r_ind,c_ind);
        if sim >= 0.9
            w=1+xita;
            M(i,3) = 1;
            W(c_ind,r_ind) = w;
            W(r_ind,c_ind) = w;
        else
            w=1-xita;
            M(i,3) = 1;
            W(c_ind,r_ind) = w;
            W(r_ind,c_ind) = w;
        end
    end
    for i = 1:pairnum_C
        c_ind = C(i,1);
        r_ind = C(i,2);
        sim = GOsim(r_ind,c_ind);
        if sim <= 0.1
            w=-1-xita;
            C(i,3) = -1;
            W(c_ind,r_ind) = w;
            W(r_ind,c_ind) = w;
        else
            w=-1+xita;
            C(i,3) = -1;
            W(c_ind,r_ind) = w;
            W(r_ind,c_ind) = w;
        end
    end
end
