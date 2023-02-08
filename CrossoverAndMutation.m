function [newpop] = CrossoverAndMutation(pop,d,alpha,pc,pm)
% Params
nm=5;
pc = 0.8;
pop_size=size(pop,1);
var_min =0;
var_max=1;
Fi=zeros(pop_size,1);
crowdis=zeros(pop_size,1);
fit = zeros(2,pop_size);
P=pop(1).P;
s = size(pop(1).Pset,2);
for i=1:pop_size
    crowdis(i)=1/pop(i).cdis;
    fit(:,i) = pop(i).f;
end

for i=1:pop_size
    ri(i)=pop(i).rank;
    Fi(i)=(alpha*(1-alpha).^(ri(i)-1));
%     Fi(i)=(alpha*(1-alpha).^(ri(i)-1))*exp(-crowdis(i));%*((1-fit(i,1))/sum(fit(:,1)));%*((1-fit(i,2))/sum(fit(:,2)));
end
Pi=cumsum(Fi)/sum(Fi);

newpop=[pop;pop];
si = 0;

    for pp=1:4:pop_size
        if rand < pm
        si = si+1;
        i1=find(rand<=Pi,1,'first');
        i2=find(rand<=Pi,1,'first');
        newpop(pp) = pop(i1);
        newpop(pp+1) = pop(i2);
        newpop(pp+2) = pop(i1);
        newpop(pp+3) = pop(i2);
        %% NDX crossover
        [newpop(pp),newpop(pp+1)] = NDXcrossover(newpop(pp),newpop(pp+1),pop(i1),pop(i2),d);
        newpop(pp+2) = newpop(pp);
        newpop(pp+3) = newpop(pp+1);
        % single point cross
        if P
            Pset1 = pop(i1).Pset;
            Pset2 = pop(i2).Pset;
            [newpop(pp+2).Pset,newpop(pp+2).Pset] = cross(Pset1,Pset2);
        end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mutation
for i=1:pop_size
    L = size(newpop(i).solution,2);
        if rand < pm 
            v = newpop(i).solution;
            u = rand(1,L);
            a = zeros(1,L);
            a(u<=0.5)=(2*u(u<=0.5) + (1-2*u(u<=0.5)).*((1-v(u<=0.5)).^(nm+1))).^(1/(nm+1))-1; 
            a(u>0.5)=1-(2*(1-u(u>0.5))+2*(u(u>0.5)-0.5).*((1-(1-v(u>0.5))).^(nm+1))).^(1/(nm+1));
            v = v + a;
            v(v>1) = 1;
            v(v<0) = 0;
            newpop(i).solution = v;
        end
        Pset = newpop(i).Pset;
        mu = rand(1,s);
        PP = setdiff(1:2*s,Pset);
        need_s = sum(mu<=pm);
        index = randperm(s);
        c = PP(index(1:need_s));
        Pset(mu<=pm) = c;
end
end


function [c1,c2]= cross(set1,set2)
s = size(set1,2);
cross_num = randi(s-1,1);
c1 = [set1(1:cross_num) set2(cross_num+1:end)];
c2 = [set2(1:cross_num) set1(cross_num+1:end)];
ind_unique_S=unique(c1);
if size(ind_unique_S,2)<s
    PP=setdiff(1:2*s,ind_unique_S);
    need_s = s - size(ind_unique_S,2);
    index = randperm(size(PP,2));
    c = PP(index(1:need_s));
    c1=[ind_unique_S c];
end

ind_unique_S=unique(c2);
if size(ind_unique_S,2)<s
    PP=setdiff(1:2*s,ind_unique_S);
    need_s = s - size(ind_unique_S,2);
    index = randperm(size(PP,2));
    c = PP(index(1:need_s));
    c2=[ind_unique_S c];
end
end