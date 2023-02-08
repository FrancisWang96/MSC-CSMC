function [U_new, beta_new] = updateU(U, V, X, W, beta, eta,model)
[k, X_n] = size(U);
D = distfcm(V, X);
switch model
    case 1
        % con
        tmp = D.^(-2);
        U_FCM= tmp./(ones(k, 1)*sum(tmp));
        U_P1 = U*W;
        U_P2 = zeros(1,X_n);
        sum_inv_di = sum(tmp,1);
        % 求uic_con的第二项
        for i = 1:X_n
            c = 1;
            sum_2_u = 0;
            for f = 1:k
                sum_2_u1 = 0;
                for j = 1:X_n
                    sum_2_u1 = sum_2_u1 + W(i,j)*U(f,j);
                end
                sum_2_u = sum_2_u + sum_2_u1*tmp(f,i);
            end
            U_P2(c,i) = sum_2_u/sum_inv_di(i);
        end
        U_P2 = ones(k,1)*U_P2;
        U_new = U_FCM + (U_P1 - U_P2).*tmp;
    case 2
        %         V = mf*data./(sum(mf,2)*ones(1,size(data,2))); %new center
    case 0
        U_FCMq = ones(k,1)*(sum(D,1)./k) - D;
        if ~isempty(U)
            U_P1 = U*W;
            %             for c = 1:k
            %                 for i = 1:X_n
            %                     T_sum = 0;
            %                     for j = 1:X_n
            %                         T_sum = T_sum + W(i,j)*U(c,j);
            %                     end
            %                     U_P1(c,i) = T_sum;
            %                 end
            %             end
            U_P2 = zeros(k,X_n);
            for c = 1:k
                for ii = 1:X_n
                    T_sum = 0;
                    for f = 1:k
                        for j = 1:X_n
                            T_sum = T_sum + W(ii,j)*U(f,j);
                        end
                    end
                    U_P2(c,ii) = T_sum/k;
                end
            end
            U_P = U_P1 - U_P2;
            % estimate beta
            if sum(U_P ~= 0) == 0
                beta_new = beta;
            else
                U_FCMq1 =  mapminmax(U_FCMq,0,1);
                U_P11 = mapminmax(U_P,0,1);
                %     beta_new = beta * (sum(sum(abs(U_FCMq1)))/sum(sum(U_FCMq1 ~= 0)))/(sum(sum(abs(U_P11)))/sum(sum(U_P11 ~= 0)));
                %     beta_new = beta * (sum(sum(abs(U_FCMq)))/sum(sum(U_FCMq ~= 0)))/(sum(sum(abs(U_P)))/sum(sum(U_P~= 0)));
                beta_new = beta;
                %     if beta_new > 2.0
                %         beta_new = 2.0;
                %     end
                % beta_new = beta;
            end
            % update U
            U_new = 1/k + (U_FCMq + U_P*beta)/eta;
        else
            U_new = 1/k + U_FCMq/eta;
        end
end
U_temp = U_new;
% 解决负值
for ii = 1:X_n
    U_new(U_new(:,ii) < 0, ii) = 0;
    U_new(U_new(:,ii) > 0, ii) = U_new(U_new(:,ii) > 0, ii)./sum(U_new(U_new(:,ii) > 0, ii));
end
beta_new = beta;
if find(all(U_new==0,2))
    U_new = mapminmax(U_temp',0,1)';
end
end