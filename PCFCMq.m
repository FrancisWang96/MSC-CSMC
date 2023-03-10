function [U, Center,obj_f,beta,i] = PCFCMq(X, U, W, beta0, eta, options)
% Input
% X(n*s)    data matrix,n sample
% W(n*n)    weight matrix
% k         numeber of clusters
% beta0      paremerter
% eta       paremerter

% Output
% U(k*n)    fuzzy memberships
% Center(k*s)     cluster centroids
if nargin ~= 5 && nargin ~= 6
    error(message("fuzzy:general:errFLT_incorrectNumInputArguments"))
end
X_n = size(X, 1);
obj_f = 0.0;
% Change the following to set default options
default_options = [	% exponent for the partition matrix U
    200;	% max. number of iteration
    1e-5;	% min. amount of improvement
    0];	% info display during iteration

if nargin == 5
    options = default_options;
else
    % If "options" is not fully specified, pad it with default values.
    if length(options) < 4
        tmp = default_options;
        tmp(1:length(options)) = options;
        options = tmp;
    end
    % If some entries of "options" are nan's, replace them with defaults.
    nan_index = find(isnan(options)==1);
    options(nan_index) = default_options(nan_index);
end

max_iter = options(1);		% Max. iteration
min_impro = options(2);		% Min. improvement
display = options(3);		% Display info or not

%% Main loop
for i = 1:max_iter
    %%
    % update the centroids
    Center = U*X./(sum(U,2)*ones(1,size(X,2)));
    % calculate Distance 
    Dist = zeros(size(Center, 1), size(X, 1));
    for k = 1:size(Center, 1)
        Dist(k, :) = sqrt(sum(((X-ones(size(X, 1), 1)*Center(k, :)).^2), 2));
    end
    %% update memberships 
    [U_P,U_FCMq] = cal_U(U, Dist, W);
    % estimate beta
    beta = cal_beta(U_P,U_FCMq,beta0);
    U_new = 1/k + (U_FCMq + U_P*beta)./eta;
    %  handle negative values
    U_temp = U_new;
    for ii = 1:X_n
        if any(U_temp(:,ii) < 0)
            U_temp(U_temp(:,ii) < 0, ii) = 0;
            U_temp(U_temp(:,ii) > 0, ii) =  U_temp(U_temp(:,ii) > 0, ii)./sum( U_temp(U_temp(:,ii) > 0, ii));
        end
    end
    col_sum = sum(U_temp);
    U_new = U_temp./col_sum(ones(k, 1), :);
    T= -sum(sum(W.*((U_new')*U_new)));
    obj_fcn(i) = sum(sum(U_new.*Dist)) + eta*0.5*sum(sum(U_new.*U_new)) + beta*T*0.5;
    %     if obj_fcn(i)<0
    %         break;
    %     end
    U = U_new;
    if display
        fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
    end
    obj_f = obj_fcn(i);
    % check termination condition
    if i > 1
        if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro
            break; 
        end
    end
end
end

function Delta = cal_centroids(X,U)
%% ???X???[s,n],U???[n,k]???
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
    %% cosin distance
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

function [U_P,U_FCMq] = cal_U(U, D, W)
    %% update memberships
    % U = 1/k + (1/eta)*(U_FCMq + beta*U_P)
    % U_P = U_P1 - U_P2/k
    [k, X_n] = size(U);
    U_FCMq = ones(k,1)*(sum(D,1)./k) - D;
    if ~isempty(U)
        U_P1 = U*W;
        U_P2 = ones(k,1)*sum(U_P1);
        U_P = U_P1 -(U_P2./k);
    end
end

function beta_new = cal_beta(U_P,U_FCMq,beta0)
    %% estimate beta
    if sum(U_P ~= 0) == 0
        beta_new = beta0;
    else
        U_FCMq1 = U_FCMq;
        U_P11 =  U_P;
        beta_new = beta0 * (sum(sum(abs(U_FCMq1)))/sum(sum(U_FCMq1 ~= 0)))/(sum(sum(abs(U_P11)))/sum(sum(U_P11 ~= 0)));
        beta_new = beta0 * (sum(sum(abs(U_FCMq)))/sum(sum(U_FCMq ~= 0)))/(sum(sum(abs(U_P)))/sum(sum(U_P~= 0)));
        %     if beta_new > 2.0
        %         beta_new = 2.0;
        %     end
    end
end




