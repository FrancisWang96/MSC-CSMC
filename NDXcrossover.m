function [object1,object2]=NDXcrossover(object1,object2,p1,p2,d)
K = size(p1.U,1);

beta = rand(1,K);%均匀分布随机数
p1.tc = chrom_to_cen(p1.solution,d);
p2.tc = chrom_to_cen(p2.solution,d);

tc1 = zeros(K,d);
tc2 = zeros(K,d);
tc1(beta<=0.5,:) = 0.5.*(p1.tc(beta<=0.5,:) + p2.tc(beta<=0.5,:) + 1.481.*(p1.tc(beta<=0.5,:) -p2.tc(beta<=0.5,:)).*abs(randn(1)));
tc1(beta>0.5,:) = 0.5.*(p1.tc(beta>0.5,:) + p2.tc(beta>0.5,:) - 1.481.*(p1.tc(beta>0.5,:) - p2.tc(beta>0.5,:)).*abs(randn(1)));

tc2(beta<=0.5,:) = 0.5.*(p1.tc(beta<=0.5,:) + p2.tc(beta<=0.5,:) - 1.481.*(p1.tc(beta<=0.5,:) -p2.tc(beta<=0.5,:))*abs(randn(1)));
tc2(beta>0.5,:) = 0.5.*(p1.tc(beta>0.5,:) + p2.tc(beta>0.5,:) + 1.481.*(p1.tc(beta>0.5,:) - p2.tc(beta>0.5,:))*abs(randn(1)));
% cf2(beta<=0.5)=0.5.*((p1.solution(beta<=0.5)+p2.solution(beta<=0.5))-1.481.*(p1.solution(beta<=0.5)-p2.solution(beta<=0.5)).*abs(ND(beta<=0.5)));
% cf2(beta>0.5)=0.5.*((p1.solution(beta>0.5)+p2.solution(beta>0.5))+1.481.*(p1.solution(beta>0.5)-p2.solution(beta>0.5)).*abs(ND(beta>0.5)));
% 截断范围
tc1(tc1>1)=1;
tc2(tc2<0)=0;

object1.solution = cen_to_chrom(tc1);
object2.solution = cen_to_chrom(tc2);

end