function [pareto,Pfit]=update(pareto,K,L0,S)
n=size(pareto,1);
[m,~]=size(pareto(1).datamatrix);
Pfit=zeros(n,2);
expo=2;
c1=0.01;
c2=0;
 %% 利用FCM迭代公式计算质心
for i=1:n 
    U=pareto(i).U;
    % updata center
    [Un,tc] = updatefcmPC(pareto(i).datamatrix,K,expo,U,S,c1,c2);
    %[Un,tc] = fcm(pop(i).datamatrix,K,expo,U);
    [~,Lab] = max(Un);
    % calculate the objective function
    [J,XB]=index_PC(pareto(i).datamatrix,Un,tc,K,S,c1,c2);
    Pfit(i,1)=J;   Pfit(i,2)=XB;    
    %%% New pop(pos,cluster,cost)
    pareto(i).clu=Lab;
    pareto(i).cost=(Pfit(i,:))';
    pareto(i).U=Un;
    pareto(i).solution(1:L0)=cen_to_chrom(tc);
end
end   

function [u,tc]=updatefcmPC(X,k,weight,u,S,c1,c2)
eps=1e-5;
[m,n]=size(X);
if ~isempty(S)
    pairS=S(:,[1 2]);
    ind_unique_S=unique(pairS);
    group_pair=FindPairRelateInS(pairS,ind_unique_S);
    labels_group=cell(1,size(group_pair,2));
    group_u0=cell(1,size(group_pair,2));
    for i=1:size(group_pair,2)
        labels_group{1,i}=unique(group_pair{1,i});
        group_u0{1,i}=zeros(length(labels_group{1,i})*k,1);%%%
    end
    ind_out_S=setdiff(1:m,ind_unique_S);
else
    ind_out_S=1:m;
end
newk=k;
ite=0;
while ite<100
    ite=ite+1;
    %% compute delta计算质心
    delta=zeros(k,n); % a row is a cluster's center
    mf = u.^weight;          
    delta = mf*X./((ones(size(X,2),1)*sum(mf'))'); 
    old_delta=delta;
    if ~isempty(num_del_cluster)
        if size(u,2)==1%如果分了1簇
            break;
        else
            u(:,num_del_cluster)=[];
            delta(num_del_cluster,:)=[];
            newk=size(u,2);
            if newk==0
                pY=ones(m,1);
                return;
            end
        end
    else
        term=0;
        for i=1:newk
            term=term+norm(delta(i,:)-old_delta(i,:))^2;
        end
        if term<eps
            break;
        end
    end
    for j=1:newk
        %dj=norm(X(ind_out_S(i),:)-delta(j,:))^2;
        dj=sum((X(ind_out_S,:)-ones(size(X(ind_out_S,:),1),1)*delta(j,:)).^2,2);
        [row,~]=find (dj<eps);%如果有距离为0的将其对应的簇的隶属度设置为1
        if size(row,1)>0
            u(ind_out_S(row),:)=0;
            u(ind_out_S(row),j)=1;
            ind_out_S_new = ind_out_S;
            ind_out_S_new(row)=[];
            dj_new= dj;
            dj_new(row)=[];
            u(j,ind_out_S_new)=1./dj_new;
        else
            u(j,ind_out_S)=1./dj;
        end
    end
    u(:,ind_out_S)=u(:,ind_out_S)./(ones(newk,1)*sum(u(:,ind_out_S)));
    u=u';
    num_del_cluster=[];
        % compute u for case 2
        if ~isempty(S)
            for i=1:size(group_pair,2)
                DforQPP=zeros(length(labels_group{1,i})*newk,length(labels_group{1,i})*newk);
                tmp=1;
                for j1=1:length(labels_group{1,i})
                    for j2=1:newk
                        DforQPP(tmp,tmp)=norm(X(labels_group{1,i}(j1),:)-delta(j2,:))^2;
                        tmp=tmp+1;
                    end
                end
                for j=1:size(group_pair{1,i},1)
                    s_thispair=GetOriginS(S,group_pair{1,i}(j,:));
                    if s_thispair>0
                        ind_thispair=find(labels_group{1,i}==group_pair{1,i}(j,1),1);
                        ind_thispair2=find(labels_group{1,i}==group_pair{1,i}(j,2),1);
                        val_thispair=c1/2*s_thispair;
                        for j2=1:newk
                            DforQPP((ind_thispair-1)*newk+j2,(ind_thispair-1)*newk+j2)=DforQPP((ind_thispair-1)*newk+j2,(ind_thispair-1)*newk+j2)+val_thispair;
                            DforQPP((ind_thispair2-1)*newk+j2,(ind_thispair2-1)*newk+j2)=DforQPP((ind_thispair2-1)*newk+j2,(ind_thispair2-1)*newk+j2)+val_thispair;
                        end
                    end
                    coef=-c1*m/size(S,1)*s_thispair; % coef of ui and uj in a pair
                    ind1=find(labels_group{1,i}==group_pair{1,i}(j,1),1);
                    ind2=find(labels_group{1,i}==group_pair{1,i}(j,2),1);
                    row_ind=(ind1-1)*newk+1:ind1*newk;
                    col_ind=(ind2-1)*newk+1:ind2*newk;
                    for r=1:length(row_ind)
                        DforQPP(row_ind(r),col_ind(r))=coef;
                        DforQPP(col_ind(r),row_ind(r))=coef;
                    end
                end
                if length(labels_group{1,i})<=3 && newk <=3
                    group_u0{1,i} = QPP_IDsolver(DforQPP,length(labels_group{1,i}),newk);
                else
                    if ~isempty(num_del_cluster)
                        group_u0{1,i}=zeros(length(labels_group{1,i})*newk,1);
                        group_u0{1,i}=QPP_IDsolver(DforQPP,length(labels_group{1,i}),newk,group_u0{1,i});
                    else
                        group_u0{1,i}=QPP_IDsolver(DforQPP,length(labels_group{1,i}),newk,group_u0{1,i});
                    end
                end
                for j=1:length(labels_group{1,i})
                    u(labels_group{1,i}(j),:)=group_u0{1,i}((j-1)*newk+1:j*newk)';
                end
            end
        end

%pY=CorrectLabel(pY);
end
tc=delta;
u=u';
end

function val = GetOriginS(S,pair)
% get the 3-rd value for the pair from S
for i=1:size(S,1)
    if S(i,1)==pair(1) && S(i,2)==pair(2)
        val=S(i,3);
        return;
    end
end
end
function [U, center] = fcm(data,cluster_n, expo,U)
mf = U.^expo;      
center = mf*data./((ones(size(data, 2), 1)*sum(mf'))');
dist = distfcm(center, data); 
tmp = dist.^(-2/(expo-1));  
U = tmp./(ones(cluster_n, 1)*sum(tmp)); 
end
function out = distfcm(center, data)
out = zeros(size(center, 1), size(data, 1));
for k = 1:size(center, 1)
    out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
end
end
function pop = cen_to_chrom(tc);
[k,d]=size(tc);
pop=[];
for i=1:k
   pop=[pop,tc(i,1:d)];
end
end
function [jm,xb]=index_PC(points,u,tc,K,S,c1,c2)
dist=distfcm(tc,points);
mf = u.^2;
sumc = sum(sum((dist.^2).*mf));
jm1=sumc;
jm2=0;
for i=1:size(S,1)
    up = S(i,1);
    uq = S(i,2);
    sval = S(i,3);
    if sval>=0
       jm2 =jm2+sum(0.5*sval*(u(:,up)-u(:,uq)).^2);
    else
       jm2 = jm2+sum(-1.0*sval*dot(u(:,up),u(:,uq)));
    end
end
jm=jm1+jm2*c1;  
disttc=distfcm(tc,tc).^2;
for i=1:K
    disttc(i,i)=9999;
end
mnm=min(min(disttc));

den=size(points,1)*mnm;
xb=sumc/den; 
% xb=1/mnm;  
kk = mean(distfcm(tc,mean(points)).^2);
kindex = (sumc+kk)/den;
end

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

