function tc=chrom_to_cen(pop,d)
[~,L]=size(pop);
Ki = L/d;
tc=zeros(Ki,d);
for i=1:Ki
   tc(i,:)= pop((i-1)*d+1:i*d);
end
end