function [spop,pareto] = Elitismselect(Rpop,NP)   
    % Initialization
    Rrank=[];
    for i=1:size(Rpop)
        Rrank=[Rrank;Rpop(i).rank];
    end
    if size(Rpop,1)<NP
        N = size(Rpop,1)
    else
        N = NP;
    end
    spop=Rpop(1:N);
    Npf = length(unique(Rrank));
    F{1}=[];
    % Selecting the chromosomes
    pf = 1;
    numberOfSolutions = 0;
    while pf <= Npf
        index=[];
        for i=1:length(Rpop)
            if Rpop(i).rank==pf
                index=[index;i];
            end
        end
        if numberOfSolutions + length(index)<= N
            spop(numberOfSolutions+1:numberOfSolutions+length(index))=Rpop(index);
            numberOfSolutions = numberOfSolutions + length(index);
        else
            rest = N - numberOfSolutions;
            temp_pop=Rpop(index);
            lastPFdist=[];
            for i=1:length(temp_pop)
                lastPFdist=[lastPFdist;temp_pop(i).cdis];
            end
            [~,idx] = sort(lastPFdist,'descend');
            temp_pop=temp_pop(idx);
            try
                spop(numberOfSolutions+1:numberOfSolutions+rest)=temp_pop(1:rest);
            catch
                display('');
            end
            numberOfSolutions = numberOfSolutions + rest;
        end
        pf = pf + 1;
    end
       
        for i=1:N
            if spop(i).rank==1
                F{1}=[F{1},i];
            end
        end
        pareto=spop(F{1},:);
end
