%GOclusterprofile(DataSetName)
function GOclusterprofile(DataSetName)
switch lower(DataSetName)
    case 'a'
        data=csvread('./Datasets/Aclus.csv',0,0);
        clusters=data(:,1);
        Data=data;
        Data(:,1)=[];
        [Y,PS] = mapminmax(Data,0,1);
        ave=mean(Y);
        va=var(Y);
        me=median(Y);
        s=std(Y);
        rng('default')
        t=0:1:7;
        figure(1)
        for c = 1:4
            subplot(2,2,c);
            plot(t,Y((clusters == c),:)','g');
            title(['cluster ',num2str(c)])
            hold on;
            ave=mean(Y((clusters == c),:));
            s=std(Y((clusters == c),:));
            plot(t,ave,'-o');
            errorbar(t, ave, s,'k');
            axis tight
        end
    case'yeast384'
        data=importdata('./Datasets/384clus.csv');
        clusters=data(:,1);
        Data=data;
        Data(:,1)=[];
        [Y,PS] = mapminmax(Data,0,1);
        ave=mean(Y);
        va=var(Y);
        me=median(Y);
        s=std(Y);
        rng('default')
        t=0:1:16;
        figure(1)
        for c = 1:5
            subplot(3,2,c);
            plot(t,Y((clusters == c),:)','g');
            title(['cluster ',num2str(c)])
            hold on;
            ave=mean(Y((clusters == c),:));
            s=std(Y((clusters == c),:));
            plot(t,ave,'-o');
            errorbar(t, ave, s,'k');
            axis tight
        end
    case'yeast237'
        data=importdata('./Datasets/237clus.csv');
        clusters=data(:,1);
        Data=data;
        Data(:,1)=[];
        [Y,PS] = mapminmax(Data,0,1);
        ave=mean(Y);
        va=var(Y);
        me=median(Y);
        s=std(Y);
        rng('default')
        t=0:1:16;
        figure(1)
        for c = 1:4
            subplot(2,2,c);
            plot(t,Y((clusters == c),:)','g');
            title(['cluster ',num2str(c)])
            hold on;
            ave=mean(Y((clusters == c),:));
            s=std(Y((clusters == c),:));
            plot(t,ave,'-o');
            errorbar(t, ave, s,'k');
            axis tight
        end
        
    case'serum'
        data=importdata('./Datasets/serumclus.csv');
        clusters=data(:,1);
        Data=data;
        Data(:,1)=[];
        [Y,PS] = mapminmax(Data,0,1);
        ave=mean(Y);
        va=var(Y);
        me=median(Y);
        s=std(Y);
        rng('default')
        t=0:1:12;
        figure(1)
        for c = 1:6
            subplot(3,2,c);
            plot(t,Y((clusters == c),:)','g');
            title(['cluster ',num2str(c)])
            hold on;
            ave=mean(Y((clusters == c),:));
            s=std(Y((clusters == c),:));
            plot(t,ave,'-o');
            errorbar(t, ave, s,'k');
            axis tight
        end
end    
end
