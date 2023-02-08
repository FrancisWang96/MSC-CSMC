function [Xdata,k,Y]=DataSelect(DataSetName,model)
Y=[];
% model:0:按特征即列对数据进行0-1标准化
% 1:对行进行z变换
% 2：
switch (DataSetName)
    case 'a'
        Xdata=importdata('./Datasets/a/rawdata/a.txt');
        Xdata = Xdata.data;
        k=4;
    case'cell'
        Xdata=importdata('./Datasets/cell/rawdata/cell.txt');
        Xdata = Xdata.data;
        Y = Xdata(:,1);
        Xdata = Xdata(:,2:18);
        k=5;
    case 'sporulation'
        Xdata =importdata('./Datasets/sporulation/rawdata/sporulation.txt');
        Xdata = Xdata.data;
        k=6;
    case'serum'
        Xdata =importdata('./Datasets/serum/rawdata/serum.txt');
        Xdata = Xdata.data;
        k=6;
    case'gal'
        Xdata =importdata('./Datasets/gal/rawdata/gal.txt');
        Xdata = Xdata.data;
        Y = Xdata(:,1);
        Xdata(:,1) = [];
        k=4;
end
if model == 0
    Xdata = mapminmax(Xdata',0,1)';
end
if model == 1
    Xdata = normalize(Xdata')';
end
if model == 2
    Xdata = normalize(Xdata')';
    Xdata = mapminmax(Xdata',0,1)';
end
end
% if model == '0'
%     switch (DataSetName)
%         case 'A'
%             Xdata=importdata('./Datasets/rawdata/A.txt');
%             Xdata = Xdata.data;
%             k=4;
%         case'yeast384'
%             Xdata=importdata('./Datasets/rawdata/yeast384.txt');
%             Xdata = Xdata.data;
%             Y = Xdata(:,1);
%             Xdata = Xdata(:,2:18);
%             Xdata = mapminmax(Xdata',0,1)';
%             k=5;
%         case 'sporulation'
%             Xdata =importdata('./Datasets/rawdata/sporulation.txt');
%             Xdata = Xdata.data;
%             Xdata = mapminmax(Xdata',0,1)';
%             k=6;
%         case'T'
%             Xdata=importdata('./Datasets/rawdata/T_dataset.txt');
%             Xdata = mapminmax(Xdata',0,1)';
%             k=6;
%         case'serum'
%             Xdata =importdata('./Datasets/rawdata/serum501.txt');
%             Xdata = Xdata.data;
%             Xdata = mapminmax(Xdata',0,1)';
%             k=6;
%         case'Gal'
%             Xdata =importdata('./Datasets/rawdata/Gal.txt');
%             Xdata = Xdata.data;
%             Y = Xdata(:,1);
%             Xdata(:,1) = [];
%             Xdata = mapminmax(Xdata',0,1)';
%             k=4;
%         case'Iris'
%             Xdata =importdata('./Datasets/rawdata/iris.txt');
%             Xdata = Xdata;
%             Y = Xdata(:,end);
%             Xdata(:,end) = [];
%             Xdata = mapminmax(Xdata',0,1)';
%             k=3;
%         case'breast'
%             Xdata =importdata('./Datasets/rawdata/breast.txt');
%             Xdata = Xdata.data;
%             Y = Xdata(:,1);
%             Xdata(:,1) = [];
%             Xdata = mapminmax(Xdata',0,1)';
%             k=6;
%     end
% end
