function  pop = Crowdingdistance(pop,F)
C=[pop.f]';
nobj=length(pop(1).f);
NF=length(F);

for i=1:NF
    NFM=length(F{i});
    C0=C(F{i},:);
    D=zeros(NFM,nobj);   
    for j=1:nobj
        Cj=C0(:,j);
        [value,index]=sort(Cj);
        minc=value(1);
        maxc=value(end);
        D(index(1),j)=Inf;
        D(index(end),j)=Inf;
        
        for k=2:NFM-1
            if maxc-minc==0
                D(index(k),j)=Inf; 
            else
                D(index(k),j)=abs(value(k+1)-value(k-1))/(maxc-minc); 
                if isnan(D(index(k),j))
                    value(k+1);
                    value(k-1);
                    maxc-minc;
                end
            end
        end
    end
    
    for z=1:NFM
       pop(F{i}(z)).cdis=sum(D(z,:)); 
    end
    
end
