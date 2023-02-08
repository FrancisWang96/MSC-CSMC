% 读取数据集
figure(1)
rng(9,'twister');
name = 'gal';
[Data,k,Y] = DataSelect(name,0);
num = 150;
[M,C,ini_P] =  pcswitch(name,num);
%画出初始的约束分布
P = ini_P;
m=size(Data,1);
A1=zeros(m,m);%ml-real
A0=zeros(m,m);%cl-real
A2=zeros(m,m);%ml-error
A3=zeros(m,m);%cl-error
for s = 1:size(P,1)
    i=P(s,1);
    j=P(s,2);
    w=P(s,4);
    if w == 1
        if P(s,3) == 1
            A1(i,j)=1;
        else
            A0(i,j)=1;
        end
    else
        if  P(s,3) == 1
            A2(i,j)=1;
        else
            A3(i,j)=1;
        end
    end
end
X=Data;
t1=31;
t2=59;
dt=[t1,t2];
X1=X(Y==1,dt);
X2=X(Y==2,dt);
X3=X(Y==3,dt);
X4=X(Y==4,dt);
scatter(X1(:,1),X1(:,2),30,'b','filled')
hold on
scatter(X2(:,1),X2(:,2),30,'c','filled')
hold on
scatter(X3(:,1),X3(:,2),30,'y','filled')
hold on
scatter(X4(:,1),X4(:,2),30,'m','filled')
box on
hold on
xlabel('维度31')
ylabel('维度59')
gplot(A0,X(:,dt),'r')
gplot(A1,X(:,dt),'g')
gplot(A2,X(:,dt),'g--')
gplot(A3,X(:,dt),'r--')

s=20;
figure(2)

c=1
subplot(2,2,c);
P = importdata("./Datasets/P/pck_p.txt");
m=size(Data,1);
A1=zeros(m,m);%ml-real
A0=zeros(m,m);%cl-real
A2=zeros(m,m);%ml-error
A3=zeros(m,m);%cl-error
sum(P(:,3) == -1)
for s = 1:size(P,1)
    i=P(s,1);
    j=P(s,2);
    w=P(s,4);
    if w == 1
        if P(s,3) == 1
            A1(i,j)=1;
        else
            A0(i,j)=1;
        end
    else
        if  P(s,3) == 1
            A2(i,j)=1;
        else
            A3(i,j)=1;
        end
    end
end



scatter(X1(:,1),X1(:,2),10,'b','filled')
hold on
scatter(X2(:,1),X2(:,2),10,'c','filled')
hold on
scatter(X3(:,1),X3(:,2),10,'y','filled')
hold on
scatter(X4(:,1),X4(:,2),10,'m','filled')
box on
hold on
xlabel('维度31')
ylabel('维度59')
gplot(A0,X(:,dt),'r')
gplot(A1,X(:,dt),'g')
gplot(A2,X(:,dt),'g--')
gplot(A3,X(:,dt),'r--')
title(['(a)'])


c=2;
subplot(2,2,c);
P = importdata("./Datasets/P/mpck_p.txt");
m=size(Data,1);
A1=zeros(m,m);%ml-real
A0=zeros(m,m);%cl-real
A2=zeros(m,m);%ml-error
A3=zeros(m,m);%cl-error
sum(P(:,3) == -1)
for s = 1:size(P,1)
    i=P(s,1);
    j=P(s,2);
    w=P(s,4);
    if w == 1
        if P(s,3) == 1
            A1(i,j)=1;
        else
            A0(i,j)=1;
        end
    else
        if  P(s,3) == 1
            A2(i,j)=1;
        else
            A3(i,j)=1;
        end
    end
end

scatter(X1(:,1),X1(:,2),10,'b','filled')
hold on
scatter(X2(:,1),X2(:,2),10,'c','filled')
hold on
scatter(X3(:,1),X3(:,2),10,'y','filled')
hold on
scatter(X4(:,1),X4(:,2),10,'m','filled')
box on
hold on
xlabel('维度31')
ylabel('维度59')
gplot(A0,X(:,dt),'r')
gplot(A1,X(:,dt),'g')
gplot(A2,X(:,dt),'g--')
gplot(A3,X(:,dt),'r--')
title(['(b)'])


c=3;
subplot(2,2,c);
P = importdata("./Datasets/P/mei_p.txt");
m=size(Data,1);
A1=zeros(m,m);%ml-real
A0=zeros(m,m);%cl-real
A2=zeros(m,m);%ml-error
A3=zeros(m,m);%cl-error
sum(P(:,3) == -1)
for s = 1:size(P,1)
    i=P(s,1);
    j=P(s,2);
    w=P(s,4);
    if w == 1
        if P(s,3) == 1
            A1(i,j)=1;
        else
            A0(i,j)=1;
        end
    else
        if  P(s,3) == 1
            A2(i,j)=1;
        else
            A3(i,j)=1;
        end
    end
end
scatter(X1(:,1),X1(:,2),10,'b','filled')
hold on
scatter(X2(:,1),X2(:,2),10,'c','filled')
hold on
scatter(X3(:,1),X3(:,2),10,'y','filled')
hold on
scatter(X4(:,1),X4(:,2),10,'m','filled')
box on
hold on
xlabel('维度31')
ylabel('维度59')
gplot(A0,X(:,dt),'r')
gplot(A1,X(:,dt),'g')
gplot(A2,X(:,dt),'g--')
gplot(A3,X(:,dt),'r--')
title(['(c)'])

s=20;
index = randperm(150);
P = ini_P(index(1:2*s),:);
X=Data;
t1=31;
t2=59;
dt=[t1,t2];
X1=X(Y==1,dt);
X2=X(Y==2,dt);
X3=X(Y==3,dt);
X4=X(Y==4,dt);

newcolors = [0 0.4470 0.7410
    %             0.8500 0.3450 0.0980
    0.50 0.16 0.16
    0.9290 0.6940 0.1250

    0.4940 0.1840 0.5560];

colororder(newcolors)

scatter(X1(:,1),X1(:,2),10,'filled')
hold on
scatter(X2(:,1),X2(:,2),10,'filled')
hold on
scatter(X3(:,1),X3(:,2),10,'filled')
hold on
scatter(X4(:,1),X4(:,2),10,'filled')
box on
hold on
xlabel('维度31')
ylabel('维度59')
%         m=size(P,1);
%         A1=zeros(m,m);
%         A0=zeros(m,m);
%         for s = 1:size(P,1)
%             i=P(s,1);
%             j=P(s,2);
%             w=P(s,4)
%             if w== 1
%                     A1(i,j)=1;
%             end
%                 if  w == -1
%                     A0(i,j)=1;
%                 end
%         end
%         gplot(A1,X(:,dt),'g--')
gplot(A0,X(:,dt),'r')
gplot(A1,X(:,dt),'g')
gplot(A2,X(:,dt),'g--')
gplot(A3,X(:,dt),'r--')


sum(P(:,4)==-1)









%
%
%
%
%
%
%
%
%
%
%         m=size(P1,1);
%         A1=zeros(m,m);
%         A0=zeros(m,m);
%
%         for s = 1:size(P1,1)
%             i=P1(s,1);
%             j=P1(s,2);
%             w=P1(s,3);
%             if w==1
%                 if label(i)==label(j)
%                     A1(i,j)=1;
%                 end
%             end
%             if w==-1
%                 if label(i)~=label(j)
%                     A0(i,j)=1;
%                 end
%             end
%         end
%         gplot(A1,X(:,dt),'g--')
%         gplot(A0,X(:,dt),'r--')
%
% %
% % for t1=1:7
% %     for t2=(t1+1):7
% %         figure()
% %         dt=[t1,t2];
% %         X1=X(class{1,1},dt);
% %         X2=X(class{1,2},dt);
% %         X3=X(class{1,3},dt);
% %         X4=X(class{1,4},dt);
% %         X5=X(class{1,5},dt);
% %         X6=X(class{1,6},dt);
% %
% %         newcolors = [0 0.4470 0.7410
% %             0.8500 0.3250 0.0980
% %             0.9290 0.6940 0.1250
% %             0.4940 0.1840 0.5560
% %             0.4660 0.6740 0.1880
% %             0.0350 0.2780 0.3540];
% %
% %         colororder(newcolors)
% %
% %
% %         scatter(X1(:,1),X1(:,2),30,'filled')
% %         hold on
% %         scatter(X2(:,1),X2(:,2),30,'filled')
% %         hold on
% %         scatter(X3(:,1),X3(:,2),30,'filled')
% %         hold on
% %         scatter(X4(:,1),X4(:,2),30,'filled')
% %         hold on
% %         scatter(X5(:,1),X5(:,2),30,'filled')
% %         hold on
% %         scatter(X6(:,1),X6(:,2),20,'filled')
% %
% %
% % scatter(X(:,t1),X(:,t2),30,'filled','bl')
% % box on
% %         hold on
% %         m=size(P,1);
% %         A1=zeros(m,m);
% %         A0=zeros(m,m);
% %         for s = 1:size(P,1)
% %             i=P(s,1);
% %             j=P(s,2);
% %             w=P(s,4)
% %             if w== 1
% %                     A1(i,j)=1;
% %             end
% %                 if  w == -1
% %                     A0(i,j)=1;
% %                 end
% %         end
% % %         gplot(A1,X(:,dt),'g--')
% %         gplot(A0,X(:,dt),'r--')
