function [outlabel,data,S] = finalcandm(pareto,K,DataSetName,percent)
[n,d]=size(pareto(1).data);
sz=length(pareto);
SPD=(0:0.05:1);
y = randsample(n,ceil(0.1*n))';
for p=1:sz
    pareto(p).project=zeros(n,d);
    for i=1:n
        for j=1:d
            for t=1:20
                if pareto(p).data(i,j)>=SPD(t) && pareto(p).data(i,j)<SPD(t+1)
                    pareto(p).project(i,j)=t;
                    continue;
                end
            end
        end
    end
    index0=pareto(p).project==0;
    pareto(p).project(index0)=20;
end
%% 
PSVIndex=zeros(sz,1);
for i=1:sz
    Lab=pareto(i).clu;
    SPDis=zeros(K,d);
    for k=1:K
        projectk=pareto(i).project(Lab==k,:);
        nk=size(projectk,1);
        for ii=1:nk
            for jj=ii+1:nk
                SPDis(k,:)=SPDis(k,:)+log(abs(projectk(ii,:)-projectk(jj,:))+1);
            end
        end
    end
    PSVIndex(i)=sum(sum(SPDis));
end
[~,kk]=min(PSVIndex);
outlabel=(pareto(kk).clu)';
data=pareto(kk).data;
S = silhouette(data,outlabel);
S=mean(S);
end

