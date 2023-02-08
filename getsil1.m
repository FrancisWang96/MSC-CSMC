function Sil = getsil(pareto,data)
[ n, ~ ] = size( data );

for i=size(pareto,1)
    [~,clust] = max(pareto(i).U);
    SI(i) = mean(silhouette(Data,clust));
end
Sil = max(SI);
end