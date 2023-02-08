function [U, V,obj_fcn,beta,i,T] = PCFCM(X, W, k, beta0, eta, options)
% Input
% X(n*s)    data matrix,n个样本s维特征
% W(n*n)    weight matrix
% k         numeber of clusters
% beta0      paremerter
% eta       paremerter
% Output
% U(k*n)    fuzzy memberships
% Delta(k*s)     cluster centroids
if nargin ~= 5 && nargin ~= 6
    error(message("fuzzy:general:errFLT_incorrectNumInputArguments"))
end
X_n = size(X, 1);
obj_f = 0.0;
% Change the following to set default options
default_options = [2;	% exponent for the partition matrix U
    200;	% max. number of iteration
    1e-5;	% min. amount of improvement
    1;	% info display during iteration
    0];

if nargin == 5
    options = default_options;
else
    % If "options" is not fully specified, pad it with default values.
    if length(options) < 5
        tmp = default_options;
        tmp(1:length(options)) = options;
        options = tmp;
    end
    % If some entries of "options" are nan's, replace them with defaults.
    nan_index = find(isnan(options) == 1);
    options(nan_index) = default_options(nan_index);
end

expo = options(1);		% Exponent for U
max_iter = options(2);		% Max. iteration
min_impro = options(3);		% Min. improvement
display = options(4);		% Display info or not
model = options(5);
%% Initialize membership  martix U
U = rand(k, X_n);
col_sum = sum(U);
U = U./col_sum(ones(k, 1), :);
beta = beta0;
%% Main loop
for i = 1:max_iter
    %         % update the centroids
    %         Delta = cal_centroids(X,U);
    %         % update memberships
    %         D = cal_D(X, Delta);% 余弦距离
    %% 欧式距离
    %% new center
    V = center(U,X,model);
    [U_new, beta_new] = updateU(U, V, X, W, beta, eta,model);
    [row,col] =  find(U_new < 0);
    Dij = pdist2(X,X,'euclidean');
    T(i) = -sum(sum(W.*((U_new')*U_new)));
    D = distfcm(V, X);
    switch model
        case 0
            obj_fcn(i) = sum(sum(U_new.*D)) + eta*0.5*sum(sum(U_new.*U_new)) + beta_new*T(i)*0.5;
        case 1
            obj_fcn(i) = sum(sum((D.^2).*(U.^2))) + beta_new*T(i);
    end

    if i > 2
        if abs(obj_fcn(i) - obj_fcn(i-1)) < 1e-3
            if abs(T(i) - T(i-1)) < 1e-3
                break;
            end
        end
    end
    beta = beta_new;
    U = U_new;
    if display
        fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
    end
    obj_f = obj_fcn(i);
    % check termination condition

    if i > 1
        if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end
    end
end
end

function Delta = cal_centroids(X,U)
%当X是[s,n],U是[n,k]时,用来验证的，现在无用
% [n,k] = size(U);
% s = size(X,1);
%Delta  =  zeros(s, k);
% for c = 1:k
%     sum_1 = 0;
%     sum_2 = 0;
%     for i = 1:n
%         sum_1 = sum_1 + U(i,c).*X(:,i);
%         sum_2 = sum_2 + U(i,c).*(X(:,i)');
%     end
%     sum_2 = sum_2*sum_1;
%     Delta(:,c) = sum_1./sqrt(sum_2);
% end
mf = U*X;
Delta = mf./(sqrt(sum(mf.*mf,2))*ones(1,size(X,2)));
end

function D_cos = cal_D(X, Delta)
% 计算数据向量与质心间的余弦距离
% D_cos(k,n)
[n,s] = size(X);
[k,~] = size(Delta);
D = zeros(k,n);
for c = 1:k
    norm_Delta_c =  Delta(c,:)/norm(Delta(c,:));
    for i = 1:n
        norm_Xi = X(i,:)/norm(X(i,:));
        D(c,i) = 1-norm_Xi*(norm_Delta_c)';
    end
end
D_cos = D;
end

function V = center(U,data,model)
expo = 2;
mf = U.^expo;       % MF matrix after exponential modification
V = [];
switch model
    case 0
        % q方式
        V = U*data./(sum(U,2)*ones(1,size(data,2))); %new center
    case 1
        % con方式
        V = mf*data./(sum(mf,2)*ones(1,size(data,2))); %new center
    case 2
        % dist方式
        V = mf*data./(sum(mf,2)*ones(1,size(data,2))); %new center
end
end






