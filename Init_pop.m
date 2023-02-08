function pop = Init_pop(X,P,pop_size,beta0,k,pareto,eta)
%随机给每个数据分标签，计算每一类数据的质心
rng(2,'twister');
[n,d]=size(X);
%% Initialize pop
empty.solution=[];%个体编码，长度为d*k
empty.f=[];
empty.dcount=[];   % dominate count
empty.dset=[];     % dominate set
empty.rank=[];
empty.cdis=inf;      % crowding distance
empty.clu=[];
empty.U=[];
empty.Evaluation = [];
empty.P = P;
s = int32(size(P,1)/2);
empty.Pset = s;
% empty.param=[];
pop=repmat(empty,pop_size,1);
var_max=1;
var_min=0;


%% Initialize centers
dist = pdist2(X, X);

%% 密度峰值聚类选出聚类中心

for i=1:round(pop_size/2)
    init_pop_i=zeros(1,k*d);
    Delta=[];
    percNeigh = 0.01+0.03*rand();
    kernel = 'Gauss';
    [dc, rho] = paraSet(dist, percNeigh, kernel);
    [Delta] = densityClust1(X,dist, dc, rho,k);
    init_pop_i=cen_to_chrom(Delta);
    [U,~] = update_FCM_U(k,Delta,X);
    %         U = calU(tc,pop(i));
    U(find(isnan(U)))=1;
    [~,Lab] = max(U);
    %     pop(i).Delta = Delta;
    pop(i).solution=init_pop_i;
    pop(i).clu=Lab;
    pop(i).U=U;
end

% 随机生成聚类中心
for i=(round(pop_size/2)+1):pop_size
    solution = [rand(d*k,1)]';
    %     center = Xsample(X,k);
    %     U = rand(k, size(X,1));
    %     col_sum = sum(U);
    %     U = U./col_sum(ones(k, 1), :);
    %     [~,Lab] = max(U);
    %     pop(i).clu=Lab;

    %     pop(i).solution=cen_to_chrom(center);
    tc = chrom_to_cen(solution,d);
    expo=2;
    dist = distfcm(tc, X);       % fill the distance matrix
    tmp = dist.^(-2/(expo-1));      % calculate new U, suppose expo != 1
    U = tmp./(ones(k,1)*sum(tmp));
    mf = U.^expo;
    center = mf*X./(sum(mf,2)*ones(1,size(X,2))); %new center
    pop(i).U = U;
    pop(i).solution   = cen_to_chrom(center);
end
% 初始化约束集
for i = 1:pop_size
    index = randperm(2*s);
    pop(i).Pset = index(1:s);
end
end

function [U,dist] = update_FCM_U(cluster_n,tc,X)
expo=2;
dist = distfcm(tc, X);
tmp = dist.^(-2/(expo-1));
U = tmp./(ones(cluster_n, 1)*sum(tmp));
end
%%
function out = distfcm(center, X)
out = zeros(size(center, 1), size(X, 1));
for k = 1:size(center, 1)
    out(k, :) = sqrt(sum(((X-ones(size(X,1),1)*center(k,:)).^2)',1));
end
end
%%
function pop=cen_to_chrom(tc)
[k,d]=size(tc);
pop=[];
for i=1:k
    pop=[pop,tc(i,1:d)];
end
end


function [centroid] = densityClust(X,dist,dc,rho,K)
[NE, ~] = size(dist);
delta = zeros(1, NE);
indNearNeigh = Inf * ones(1, NE);
[~, ordRho] = sort(rho, 'descend');

for i = 2 : NE
    delta(ordRho(i)) = max(dist(ordRho(i), :));
    for j = 1 : (i-1)
        if dist(ordRho(i), ordRho(j)) < delta(ordRho(i))
            delta(ordRho(i)) = dist(ordRho(i), ordRho(j));
            indNearNeigh(ordRho(i)) = ordRho(j);
        end
    end
end
delta(ordRho(1)) = max(delta);
indNearNeigh(ordRho(1)) = 0;
for i=1:NE
    ind(i)=i;
    gamma(i)=rho(i)*delta(i);
end
[val, Index] = sort(gamma, 'descend');
if val(K)-val(K+1)<0.2
    if rand(1)<0.7
        icl(1:K-1)=Index(1:K-1);
        icl(K)=Index(K+1);
    end
else
    icl = Index(1:K);
end
icl = Index(1:K);
centroid=[];
for i=1:K
    centroid=[centroid;X(icl(i),:)];
end
end

function [centroid] = densityClust1(X,dist,dc,rho,K)
[NE, ~] = size(dist);
delta = zeros(1, NE);
indNearNeigh = Inf * ones(1, NE);
[~, ordRho] = sort(rho, 'descend');

for i = 2 : NE
    delta(ordRho(i)) = max(dist(ordRho(i), :));
    for j = 1 : (i-1)
        if dist(ordRho(i), ordRho(j)) < delta(ordRho(i))
            delta(ordRho(i)) = dist(ordRho(i), ordRho(j));
            indNearNeigh(ordRho(i)) = ordRho(j);
        end
    end
end
delta(ordRho(1)) = max(delta);
indNearNeigh(ordRho(1)) = 0;
for i=1:NE
    ind(i)=i;
    gamma(i)=rho(i)*delta(i);
end
[val, Index] = sort(gamma, 'descend');
%   icl = Index(1:K);
%if val(K)-val(K+1)<0.2
piece = randperm(K+3);
ic=Index(piece(1:K));
%icl(1icl:K-1)=Index(1:K-1);
%icl(K)=Index(K+1);
% end

%icl = Index(1:K);
centroid=[];
for i=1:K
    centroid=[centroid;X(ic(i),:)];
end
end

function [dc, rho] = paraSet(dist, percNeigh, kernel)
distRow = squareform(dist, 'tovector');
sortDistRow = sort(distRow);
[NE, ~] = size(dist);
dc = sortDistRow(round((NE*(NE-1)/2)*percNeigh));

if strcmp(kernel, 'Gauss')
    rho = sum(exp(-(dist/dc).^2) , 1) - exp(-(0/dc).^2);
elseif strcmp(kernel, 'Cut-off')
    rho = sum(dist < dc, 1) - 1;
else

end

end