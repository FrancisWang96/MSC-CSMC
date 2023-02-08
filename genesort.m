%  BP= csvread("./Datasets/Sim_Y2genedataBP.csv",1,1);
%  Xdata=importdata('./Datasets/yeast237.txt');
%   BP= csvread("./Datasets/Sim_SgenedataBP.csv",1,1);
%  Xdata=importdata('./Datasets/serum.txt');
gene_sim =1.0- csvread("./Datasets/yeast_wang.csv",1,1);
n=size(gene_sim,1);
count = sum(gene_sim==0,2);
select_ss = gene_sim(count~=n-1,count~=n-1) ;
local = (1:n);
ind = local(count~=n-1);
select_n = size(select_ss,2);
X= csvread("./Datasets/yeast_ex.csv",1,1);
X = mapminmax(X',0,1)';
D = pdist2(X,X);
select_D = D(count~=n-1,count~=n-1);
ds = (select_D.^2)./select_ss;

nn =select_n*(select_n-1)/2;
Gij=zeros(nn,3); 

kk=0;
for i=1:select_n-1
    for j =i+1:select_n
        kk=kk+1;
        Gij(kk,1)=i;
        Gij(kk,2)=j;
        Gij(kk,3)=ds(i,j);
    end
end
% id =  isinf(Gij(:,3)); 
% [inr,inc]=find(id ==1);
% Gij(inr,:)=[];
[sortG,index] = sort(Gij(:,3));
sort_pairS=Gij(index,:);



 BP= csvread("./Datasets/Sim_AgenedataBP.csv",1);
 Xdata=importdata('./Datasets/A.txt');
 %Xdata=Xdata(:,2:13);


% dlmwrite('./Datasets/yeast237_Gij.txt',p,'-append');  
%  X=importdata('./Datasets/yeast237_Gij.txt');
dlmwrite('./Datasets/yeast_pairS.txt',sort_pairS,'-append');  
 X=importdata('./Datasets/yeast_pairS.txt');