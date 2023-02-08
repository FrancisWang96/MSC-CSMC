function U = calU(tc,ipop,data)
%% 得到当前个体的信息
S=ipop.pairS;
c1=ipop.weight(1);
c2=ipop.weight(2);
[m,n]=size(data);
k=size(tc,1);
u=zeros(m,k);
delta = tc;
%% 有成对约束时
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
%% 计算隶属度u
for j=1:k
    dj=sum((data(ind_out_S,:)-ones(size(data(ind_out_S,:),1),1)*delta(j,:)).^2,2);
    [row,~]=find (dj<eps);%如果有距离为0的将其对应的簇的隶属度设置为1
    if size(row,1)>0
        u(ind_out_S(row),:)=0;
        u(ind_out_S(row),j)=1;
        ind_out_S_new = ind_out_S;
        ind_out_S_new(row)=[];
        dj_new= dj;
        dj_new(row)=[];
        u(ind_out_S_new,j)=1./dj_new;
    else
        u(ind_out_S,j)=1./dj;
    end
end
% u=u';
% u(:,ind_out_S)=u(:,ind_out_S)./(ones(k,1)*sum(u(:,ind_out_S)));
u(ind_out_S,:)=u(ind_out_S,:)./((ones(k,1)*sum(u(ind_out_S,:)'))');
% compute u for case 2 有约束的样本的隶属度计算
if ~isempty(S)
    for i=1:size(group_pair,2)
        DforQPP=zeros(length(labels_group{1,i})*k,length(labels_group{1,i})*k);
        tmp=1;
        for j1=1:length(labels_group{1,i})
            for j2=1:k
                DforQPP(tmp,tmp)=norm(data(labels_group{1,i}(j1),:)-delta(j2,:))^2;
                tmp=tmp+1;
            end
        end
        for j=1:size(group_pair{1,i},1)
            s_thispair=GetOriginS(S,group_pair{1,i}(j,:));
            if s_thispair>0
                ind_thispair=find(labels_group{1,i}==group_pair{1,i}(j,1),1);
                ind_thispair2=find(labels_group{1,i}==group_pair{1,i}(j,2),1);
                val_thispair=c1/2*s_thispair;
                for j2=1:k
                    DforQPP((ind_thispair-1)*k+j2,(ind_thispair-1)*k+j2)=DforQPP((ind_thispair-1)*k+j2,(ind_thispair-1)*k+j2)+val_thispair;
                    DforQPP((ind_thispair2-1)*k+j2,(ind_thispair2-1)*k+j2)=DforQPP((ind_thispair2-1)*k+j2,(ind_thispair2-1)*k+j2)+val_thispair;
                end
            end
            coef=-c1*m/size(S,1)*s_thispair; % coef of ui and uj in a pair
            ind1=find(labels_group{1,i}==group_pair{1,i}(j,1),1);
            ind2=find(labels_group{1,i}==group_pair{1,i}(j,2),1);
            row_ind=(ind1-1)*k+1:ind1*k;
            col_ind=(ind2-1)*k+1:ind2*k;
            for r=1:length(row_ind)
                DforQPP(row_ind(r),col_ind(r))=coef;
                DforQPP(col_ind(r),row_ind(r))=coef;
            end
        end
        if length(labels_group{1,i})<=3 && k <=3
            group_u0{1,i} = QPP_IDsolver(DforQPP,length(labels_group{1,i}),k);
        else
%             if ~isempty(num_del_cluster)
%                 group_u0{1,i}=zeros(length(labels_group{1,i})*k,1);
%                 group_u0{1,i}=QPP_IDsolver(DforQPP,length(labels_group{1,i}),k,group_u0{1,i});
%             else
                 group_u0{1,i}=QPP_IDsolver(DforQPP,length(labels_group{1,i}),k,group_u0{1,i});
%             end
         end
        for j=1:length(labels_group{1,i})
            u(labels_group{1,i}(j),:)=group_u0{1,i}((j-1)*k+1:j*k)';
        end
    end
end
U=u';
end


function val=GetOriginS(S,pair)
% get the 3-rd value for the pair from S
for i=1:size(S,1)
    if S(i,1)==pair(1) && S(i,2)==pair(2)
        val=S(i,3);
        return;
    end
end
end