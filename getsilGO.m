function Sil_list = getsilGO(pareto_list,Data)
rng(2,"philox");
[ n, ~ ] = size( Data );
Sil_list = [];
for i = 1:size(pareo_list,2)
    pareto = pareto_list(i);
    for t=1:size(pareto,1)
        [~,clust] = max(pareto(t).U);
        SI(t) = mean(silhouette(Data,clust));
    end
    Sil = max(SI);
end
Sil_list(i) = Sil;
end
