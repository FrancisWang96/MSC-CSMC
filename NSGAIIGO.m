function Sil_list=NSGAIIGO(DataSetName,pareto,xita)
%% 种群参数设置
% alpha 轮盘赌选择参数 (alpha*(1-alpha).^(ri(i)-1))
% pm 种群变异参数
% maxgen, pop_size 种群迭代最大数个体数
model = 0;
rng(2,'twister');
% 读取数据集
[Data,k,Y] = DataSelect(DataSetName,model);
GOsim = importGO(DataSetName);
Param =[0.3 0.1 300 100 0.1 0.001];
alpha = Param(1);
pc = 0.8;
pm = Param(2);
maxgen  = Param(3);
pop_size = Param(4);
beta0 = Param(5);
eta = Param(6);
[n,d] = size(Data);
s_list=[5 10 15 20 25];
[M,C,ini_P] =  pcswitch(DataSetName,150);
pareto_list = {};
for i = 1:size(s_list,2)
    s=s_list(i);
    index = randperm(2*s);
    P = ini_P(index,:);
    gen = 1
    pop= Init_pop(Data,P,pop_size,beta0, k,pareto,eta);
    [pop,~] = FitnessGO(Data,pop,beta0,eta,GOsim,xita);
    [pop,F] = NDsort(pop);
    pop = Crowdingdistance(pop,F);
    %非支配解集pareto
    pareto = pop(find([pop.rank] == 1));

    %% 优化
    stopCondition = false;
    while ~stopCondition
        %选择,交叉，突变产生子代pop_M
        [pop_M] = CrossoverAndMutation(pop,d,alpha,pc,pm);
        Rpop=[pop;pop_M];
        [Rpop,~]=FitnessGO(Data,Rpop,beta0,eta,GOsim,xita);
        %     Rpop=[pop;pop_M];
        Rpop = check(Rpop);
        % 非支配排序
        [Rpop,F]=NDsort(Rpop);
        Rpop=Crowdingdistance(Rpop,F);
        [pop,pareto] = Elitismselect(Rpop,pop_size);
        gen = gen + 1
        % 迭代后期beta的值增大
        if(floor(maxgen*0.5) <= gen && gen <(floor(maxgen*0.5)+1))
            beta0 = beta0*2;
        end
        if(floor(maxgen*0.9) <= gen && gen <(floor(maxgen*0.9)+1))
            beta0 = beta0*2;
        end
        if(gen>=maxgen)
            stopCondition = true;
        end
    end
    pareto_list(i) = {pareto};
end
Sil_list = getsilGO(pareto_list,Data);
end

function pop=check(pop)
pop_size = size(pop,1);
index = [];
for i = 1: pop_size
    if isnan(pop(i).U)
        index=[index i];
    end
end
pop(index)=[];
end


function [sort_P] = dele_S(pareto,P)
M = P(P(:,3)>0,:);
C = P(P(:,3)<0,:);
%% 计算pareto解中满足约束的个数
for i=1:size(M,1)
    i1 = P(i,1);
    i2 = P(i,2);
    correct_num = 0;
    for j = 1:size(pareto,1)
        i1_label = pareto(j).clu(i1);
        i2_label = pareto(j).clu(i2);
        if i1_label == i2_label
            flag = 1;
        else
            flag = 0;
        end
        correct_num = correct_num +flag;
    end
    P(i,5) = correct_num;
    P(i,6) = size(pareto,1);
end
for i=(size(M,1)+1):size(P,1)
    i1 = P(i,1);
    i2 = P(i,2);
    correct_num = 0;
    for j = 1:size(pareto,1)
        i1_label = pareto(j).clu(i1);
        i2_label = pareto(j).clu(i2);
        if i1_label == i2_label
            flag = 0;
        else
            flag = 1;
        end
        correct_num = correct_num +flag;
    end
    P(i,5) = correct_num;
    P(i,6) = size(pareto,1);
end
% for i=1:size(M,1)
%     if M(i,3) > median(M(:,3))
%         if (P(i,5)) < size(pareto,1)
%             dele = [dele;i]
%         end
%     end
% end
% for i=(size(M,1)+1):size(P,1)
%     if P(i,3) < median(P((size(M,1)+1):size(P,1),3))
%         if P(i,5) < size(pareto,1)
%             dele = [dele;i]
%         end
%     end
%     P(dele,:)=[];
%     %         if isempty(dele)
%     %             P = P(1:s_list(t),:)
%     %         end
% end
[~,inds] = sort(P(:,5),'descend');
sort_P = P(inds,:);
%         if isempty(dele)
%             P = P(1:s_list(t),:)
%         end
end