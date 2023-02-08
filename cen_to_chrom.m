function pop=cen_to_chrom(tc)
[k,d]=size(tc);
pop=[];
for i=1:k
   pop=[pop,tc(i,1:d)];
end
end