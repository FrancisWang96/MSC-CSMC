function [U,V_new,beta_new] = updateU(X,V,UO, W, k,beta)
for i = 1:size(X,1)
    dic2 = sum((ones(k, 1)*X(i,:)-V).^2,2);
    Ufcm=[];
    Ufcm = (1./dic2)./sum(1./dic2,1);
    temp=0;
    for c =1:k
        temp=temp+sum(W(i,:).*UO(c,:),2)./dic2(c,:);
    end
    U(:,i) = Ufcm;
    if sum(sum(W))
        Upc= [];
        [I,y]=find(W~=0);
        if find(I==i)
            for c =1:k
                Upc(c) = (0.5*beta/dic2(c))*(sum(W(i,:).*UO(c,:),2) - temp/sum(1./dic2,1));
            end
            U(:,i) = Ufcm+Upc';
        end
    end
end
% 解决负值
for ii = 1:size(X,1)
    U(U(:,ii) < 0, ii) = 0;
    U(U(:,ii) > 0, ii) = U(U(:,ii) > 0, ii)./sum(U(U(:,ii) > 0, ii));
end
beta_new = beta;
% 得到更新后的U,更新V
mf = U.^2;       % MF matrix after exponential modification
V_new = mf*X./(sum(mf,2)*ones(1,size(X,2))); %new center
end