function [pop,Pfit]=FitnessGO(Data,pop,beta0,eta,GOsim,xita)
    pop_size = size(pop,1);
    [n,d]=size(Data);
    Pfit=zeros(pop_size,2);
    P=pop(1).P;
    %% 利用迭代公式计算质心
    for i=1:pop_size
        U=pop(i).U;
        Center = chrom_to_cen(pop(i).solution,d);
        k=size(U,1);
        Pset = pop(i).Pset;
        PairS = P(Pset,:);
        if size(Pset)
        [~,~,W] = cal_W_GO(PairS(1:3,:),Data,GOsim,xita);
        else
            W =zeros(n,n);
        end
        [f1,f2,U_new,Center_new] = fit12(Data,U,Center,W,beta0,eta,1,1e-4);
        pop(i).U = U_new;
        [~,Lab] = max(U);
        Pfit(i,1)=f1;
        Pfit(i,2)=f2;
        pop(i).clu=Lab;
        pop(i).f=(Pfit(i,:))';
        pop(i).solution(1:k*d)=cen_to_chrom(Center_new);
    end
end

function pop = cen_to_chrom(tc);
[k,d]=size(tc);
pop=[];
for i=1:k
    pop=[pop,tc(i,1:d)];
end
end


function [Jpcmq,XB,U_new,Center_new] = fit12(Data,U,Center,W,beta0,eta,max_iter,min_impro)
%% calculate objective function Jpcmq
Dist = distfcm(Center, Data);
[k,X_n] = size(U);
for i = 1:max_iter
%% update memberships 
    [U_P,U_FCMq] = cal_U(U, Dist, W);
     % estimate beta
    beta = cal_beta(U_P,U_FCMq,beta0);
    U_new = 1/k + (U_FCMq + U_P*beta)./eta;
    %  handle negative values
    U_temp = U_new;
    for ii = 1:X_n
        if any(U_temp(:,ii) < 0)
            U_temp(U_temp(:,ii) < 0, ii) = 0;
            U_temp(U_temp(:,ii) > 0, ii) =  U_temp(U_temp(:,ii) > 0, ii)./sum( U_temp(U_temp(:,ii) > 0, ii));
        end
    end
    col_sum = sum(U_temp);
    U_new = U_temp./col_sum(ones(k, 1), :);
    T= -sum(sum(W.*((U_new')*U_new)));
    obj_fcn(i) = sum(sum(U_new.*Dist)) + eta*0.5*sum(sum(U_new.*U_new)) + beta*T*0.5;
    %     if obj_fcn(i)<0
    %         break;
    %     end
    Center_new = U_new*Data./(sum(U_new,2)*ones(1,size(Data,2)));
    Dist = distfcm(Center_new, Data);
    obj_f = obj_fcn(i);
    % check termination condition
    if i > 1
        if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro
            break; 
        end
    end
end
Jpcmq = obj_f;
%% calculate objective function XB
expo=2;
% Dist = distfcm(Center, Data);
J = sum(sum((Dist.^2).*(U_new.^expo)));
% calculate D1
distc1 = distfcm(Center_new, Center_new);
for i=1:k
    distc1(i,i)=999999;
end
D1=min(min(distc1));
D1=D1.^2;
denominator=X_n*D1;
XB=J/denominator;
end


function Delta = cal_centroids(X,U)
%% 当X是[s,n],U是[n,k]时
% [n,k] = size(U);
% s = size(X,1);
%Delta  =  zeros(s, k);
% for c = 1:k
%     sum_1 = 0;
%     sum_2 = 0;
%     for i = 1:n
%         sum_1 = sum_1 + U(i,c).*X(:,i);
%         sum_2 = sum_2 + U(i,c).*(X(:,i)');
%     end
%     sum_2 = sum_2*sum_1;
%     Delta(:,c) = sum_1./sqrt(sum_2);
% end
mf = U*X;
Delta = mf./(sqrt(sum(mf.*mf,2))*ones(1,size(X,2)));
end

function D_cos = cal_D(X, Delta)
    %% cosin distance
    % D_cos(k,n)
    [n,s] = size(X);
    [k,~] = size(Delta);
    D = zeros(k,n);
    for c = 1:k
        norm_Delta_c =  Delta(c,:)/norm(Delta(c,:));
        for i = 1:n
            norm_Xi = X(i,:)/norm(X(i,:));
            D(c,i) = 1-norm_Xi*(norm_Delta_c)';
        end
    end
    D_cos = D;
end

function [U_P,U_FCMq] = cal_U(U, D, W)
    %% update memberships
    % U = 1/k + (1/eta)*(U_FCMq + beta*U_P)
    % U_P = U_P1 - U_P2/k
    [k, X_n] = size(U);
    U_FCMq = ones(k,1)*(sum(D,1)./k) - D;
    if ~isempty(U)
        U_P1 = U*W;
        U_P2 = ones(k,1)*sum(U_P1);
        U_P = U_P1 - (U_P2./k);
    end
end

function beta_new = cal_beta(U_P,U_FCMq,beta0)
    %% estimate beta
    if sum(U_P ~= 0) == 0
        beta_new = beta0;
    else
        U_FCMq1 = U_FCMq;
        U_P11 =  U_P;
        beta_new = beta0 * (sum(sum(abs(U_FCMq1)))/sum(sum(U_FCMq1 ~= 0)))/(sum(sum(abs(U_P11)))/sum(sum(U_P11 ~= 0)));
        beta_new = beta0 * (sum(sum(abs(U_FCMq)))/sum(sum(U_FCMq ~= 0)))/(sum(sum(abs(U_P)))/sum(sum(U_P~= 0)));
        %     if beta_new > 2.0
        %         beta_new = 2.0;
        %     end
    end
end


