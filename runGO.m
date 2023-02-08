rng(2,'twister');
% s = rng;
% rng(s);

 DataSetName="a";
 xita = 0.4;

% DataSetName = "gal";
% xita = 0.4;

%DataSetName="cell";
%xita = 0.6;

%DataSetName = "sporulation";
% xita = 0.7;

% DataSetName = "serum";
% xita = 0.5;
%%
pareto = [];
sil_list = NSGAIIGO(DataSetName, pareto,xita);