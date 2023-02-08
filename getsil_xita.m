function Sil= getsil_xita(pareto,Data,xita)
rng(2,"philox");

for t=1:size(pareto,1)
    [~,clust] = max(pareto(t).U);
    SI(t) = mean(silhouette(Data,clust));
end
Sil = max(SI);
end


