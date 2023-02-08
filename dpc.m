dist = csvread('H:/cell_sim_AIC.csv',1,1); 
dist = round(dist,5);
dist = 1-dist;
percNeigh = 0.02;
kernel = 'Gauss';
K=5;
[dc, rho] = paraSet(dist, percNeigh, kernel); 
[tc] = densityClust(dist, dc, rho,K);
init_pop(i,1:L0)=cen_to_chrom(tc);
[U,~] = initfcm(K,tc,data);
U(find(isnan(U)))=1;
[~,Lab] = max(U);
pop(i).solution=init_pop(i,:);
pop(i).clu=Lab;
pop(i).U=U;
%pop(i).fillsol=fillpop(i).fillsol;
%pop(i).evaluerange=fillpop(i).evaluerange;
pop(i).datamatrix=data;

for i=2:pop_size
    tc=[];
   for j=1:K
      r=round(unifrnd(1,n));
      tc(j,:)=data(r,:);
   end
   init_pop(i,1:L0)=cen_to_chrom(tc);
%     percNeigh = 0.01+0.03*rand();
%     kernel = 'Gauss';
%     [dc, rho] = paraSet(dist, percNeigh, kernel); 
%     [tc] = densityClust1(data,dist, dc, rho,K);
%     init_pop(i,1:L0)=cen_to_chrom(tc);
    [U,~] = initfcm(K,tc,data);
    U(find(isnan(U)))=1;
    [~,Lab] = max(U);
    pop(i).solution=init_pop(i,:);
    pop(i).clu=Lab;
    pop(i).U=U;
    %pop(i).fillsol=fillpop(i).fillsol;
    %pop(i).evaluerange=fillpop(i).evaluerange;
    pop(i).datamatrix=data;
end

function [U,dist] = initfcm(cluster_n,tc,data)
expo=2;
dist = distfcm(tc, data); 
tmp = dist.^(-2/(expo-1));  
U = tmp./(ones(cluster_n, 1)*sum(tmp)); 
end
%%
function out = distfcm(center, data)
out = zeros(size(center, 1), size(data, 1));
for k = 1:size(center, 1)
    out(k, :) = sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
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

function [ic] = densityClust(dist,dc,rho,K)    
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
            ic(1:K-1)=Index(1:K-1);
            ic(K)=Index(K+1);
        end
    else
        ic = Index(1:K);
    end
    ic = Index(1:K);
end

function [centroid] = densityClust1(data,dist,dc,rho,K)    
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
        centroid=[centroid;data(ic(i),:)];
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