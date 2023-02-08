label = pareto(1).clu;
class={};
for c=1:k
    class(c)={find(label==c)};
end
X1=X(class{1,1},1,3);
X2=X(class{1,2},1,3);
X3=X(class{1,3},1,3);
X4=X(class{1,4},1,3);
X5=X(class{1,5},1,3);
X6=X(class{1,6},1,3);

figure(1)
newcolors = [0 0.4470 0.7410
            0.8500 0.3250 0.0980
             0.9290 0.6940 0.1250
             0.4940 0.1840 0.5560
             0.4660 0.6740 0.1880
             0.3010 0.7450 0.9330];
         
colororder(newcolors)


scatter(X1(:,1),X1(:,2),20,'filled')
hold on
scatter(X2(:,1),X2(:,2),20,'filled')
hold on
scatter(X3(:,1),X3(:,2),20,'filled')
hold on
scatter(X4(:,1),X4(:,2),20,'filled')
hold on
scatter(X5(:,1),X5(:,2),20,'filled')
hold on
scatter(X6(:,1),X6(:,2),20,'filled')
m=size(P1,1);
A1=zeros(m,m);
A0=zeros(m,m);
for s = 1:size(P1,1)
    i=P1(s,1);
    j=P1(s,2);
    w=P1(s,3);
    if w==1
        if label(i)==label(j)
            A1(i,j)=1;
        else
            A0(i,j)=1;
        end
    end
    if w==-1
        if label(i)~=label(j)
             A1(i,j)=1;
        else
             A0(i,j)=1;
        end
    end
end
gplot(A1,X(:,1,3),'--g')
gplot(A0,X(:,1,3),'--r')