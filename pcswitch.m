function [M,C,P] = pcswitch(DataSetName,num)
M=[];
C=[];
if DataSetName =="gal"
    Xdata = importdata('./Datasets/gal/pcdata/Gal_constraints_150.txt');
else
    switch (DataSetName)
        case 'a'
            path = ['./Datasets/a/pcdata/a_constraints_',num2str(num)];
            file = [path,'.csv'];
            Xdata=csvread(file);
        case'cell'
            path = ['./Datasets/cell/pcdata/cell_constraints_',num2str(num)];
            file = [path,'.csv'];
            Xdata=csvread(file);
        case 'sporulation'
            path = ['./Datasets/sporulation/pcdata/sporulation_constraints_',num2str(num)];
            file = [path,'.csv'];
            Xdata=csvread(file);
        case'serum'
            path = ['./Datasets/serum/pcdata/serum_constraints_',num2str(num)];
            file = [path,'.csv'];
            Xdata=csvread(file);
    end

end
% X = Xdata(find(Xdata(:,3)==1),:);%选择正确约束
X = Xdata;
X(:,1:2) = X(:,1:2)+1;
M=X(find(X(:,4)==0),:);
M=[M;X(find(X(:,4)==2),:)];
C=[X(find(X(:,4)==1),:)];
M(:,4) = M(:,3);
C(:,4) = C(:,3);
a=1;
M(:,3) = a;
C(:,3) = -a;
P = [M;C];
end
