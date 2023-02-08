function [Xdata,K]=importdata(DataSetName)
switch lower(DataSetName)
    case 'a'
        Xdata=importdata('./Datasets/A.txt');
        K=4;
    case'yeast384'
%          Xdata=importdata('./Datasets/yeast384.txt');
%          Xdata = Xdata.data;
%          Xdata = Xdata(:,2:18);
         Xdata= csvread("./Datasets/yeast_ex.csv",1,1);
%          completedata = log2(completedata);
%          completedata =zscore(completedata')';

         S= csvread("./Datasets/yeast384_Gij.csv");
% 
% S(:,3)=-mapminmax(S(:,3)',-1,1)';
% SP = old_S(old_S(:,3)>0,:);%正约束
% SN = old_S(old_S(:,3)<=0,:);%负约束
% SP=SP(1:30,:);
         Xdata= mapminmax(Xdata',0,1)';
%          X = csvread('H:/CELL_EX.csv',1,1); 
%          X= mapminmax(X',0,1)';
%          completedata= X;
         K=5;
    case 'sporulation'
        Xdata =importdata('./Datasets/sporulation.txt');
        Xdata = Xdata.data;
        Xdata = mapminmax(Xdata',0,1)';
        K=6;
    case'yeast237'
         Xdata=importdata('./Datasets/yeast237.txt');
         K=4;
    case'serum'
        Xdata=importdata('./Datasets/serum.txt');
        K=6;
end
end