rng(2,'twister')
DataSetName_list = ["gal" "cell" "sporulation" "serum" "a"] ;
xita_list = 0.1:0.1:0.9;
Sil_list =zeros(5,9);
%%
for dataname = 1:5
    DataSetName=DataSetName_list(dataname);
    pareto = [];
    model = 0;
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
    Sil=[];
    for i = 1:size(xita_list,2)
        s=15;
        xita = xita_list(i);
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
        Sil_list(dataname,i) = getsil_xita(pareto,Data,xita);
    end
end

figure()
plot(xita_list,Sil_list(1,:),'b-x')
hold on
plot(xita_list,Sil_list(2,:),'r-o')
hold on
plot(xita_list,Sil_list(3,:),'g-d')
hold on
plot(xita_list,Sil_list(4,:),'ms-')
hold on
plot(xita_list,Sil_list(5,:),'k-v')
hold on
ylim([0.3,0.8]);
xlabel('\theta');
ylabel('SI');
% xlabel()
