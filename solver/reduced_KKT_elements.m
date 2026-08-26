function [ru_0,rse_k,ru_k,rse_ter,R_0,Q_k,R_k,Y_k,Q_ter,...
          g2_0,g2_k,g2_ter,v2_0,v2_k,v2_ter,...
          rv_v2_0,rv_v2_k,rv_v2_ter,ri_0,ri_k,ri_ter,...
          iS_0,iS_k,iS_ter,...
          iS_ri_hat_0,iS_ri_hat_k,iS_ri_hat_ter] = reduced_KKT_elements(mpc,...
                        t,N,nu,nx,R_0,Q_k,R_k,Y_k,Q_ter,...
                        R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Q_f0_ter,...
                        ru_0,rse_k,ru_k,rse_ter,...
                        grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
                        u_cnstr,du_cnstr,s_cnstr,s_col,su_col,ng_k,nv_k,...
                        g_0,g_k,g_ter,v_0,v_k,v_ter,g2_0,g2_k,g2_ter,v2_0,v2_k,v2_ter,...
                        rv_v2_0,rv_v2_k,rv_v2_ter,ri_0,ri_k,ri_ter,...
                        grad_qv_0,grad_qv_k,grad_qv_ter,v_rows_0,v_rows_k,...
                        iS_0,iS_k,iS_ter,iS_ri_hat_0,iS_ri_hat_k,iS_ri_hat_ter)

% f0 terms
ru_0 = grad_u_f0_0;
rse_k = grad_se_f0_k;
ru_k = grad_u_f0_k;
rse_ter = grad_se_f0_ter;

R_0 = t * R_f0_0;
Q_k = t * Q_f0_k;
R_k = t * R_f0_k;
Y_k = t * Y_f0_k;
Q_ter = t * Q_f0_ter;
%% ri_hat = ri - g^2*(-1/g) + v^2*(t*qv-1/v)
%  ri_hat = ri + g + t*qv*v^2 - v

if ng_k(1)
    ri_0 = ri_0+g_0;
    g2_0 = g_0.^2;
    if nv_k(1)
        v2_0 = v_0.^2;
        % v^2*rv = t*qv*v^2 - v
        rv_v2_0 = t*grad_qv_0.*v2_0-v_0;
        ri_0(v_rows_0) = ri_0(v_rows_0) + rv_v2_0;
    end
end

if ng_k(2)
    ri_k = ri_k+g_k;
    g2_k = g_k.^2;
    if nv_k(2)
        v2_k = v_k.^2;
        % v^2*rv = t*qv*v^2 - v
        rv_v2_k = t*grad_qv_k.*v2_k-v_k;
        ri_k(v_rows_k,:) = ri_k(v_rows_k,:) + rv_v2_k;
    end
end

if ng_k(3)
    g2_ter = g_ter.^2;
    v2_ter = v_ter.^2;
    % v^2*rv = t*qv*v^2 - v
    rv_v2_ter = t*grad_qv_ter.*v2_ter-v_ter;
    ri_ter = ri_ter+g_ter + rv_v2_ter;
end

%% S = g^2 + v^2

if ng_k(1)
    iS_0 = g2_0;
    if nv_k(1), iS_0(v_rows_0) = iS_0(v_rows_0) + v2_0; end
    iS_0  = 1./iS_0;
    iS_ri_hat_0 = iS_0.*ri_0;
end

if ng_k(2)
    iS_k = g2_k;
    if nv_k(2), iS_k(v_rows_k,:) = iS_k(v_rows_k,:) + v2_k; end
    iS_k  = 1./iS_k;
    iS_ri_hat_k = iS_k.*ri_k;
end

if ng_k(3)
    iS_ter = g2_ter;
    if nv_k(3), iS_ter = iS_ter + v2_ter; end
    iS_ter  = 1./iS_ter;
    iS_ri_hat_ter = iS_ter.*ri_ter;
end

%% rx_hat = t*grad_x + Ai*S^-1*ri_hat

if mpc.has_u_cnstr
    if u_cnstr.min_limit && u_cnstr.max_limit

        row_min = u_cnstr.min_ineqRow_0;
        row_max = u_cnstr.max_ineqRow_0;

        % k = 0
        for i = 1:nu
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            iS_min_i = iS_0(row_min_i);
            iS_max_i = iS_0(row_max_i);
            iS_ri_hat_min_i = iS_ri_hat_0(row_min_i);
            iS_ri_hat_max_i = iS_ri_hat_0(row_max_i);

            ru_0(i) = ru_0(i) + ...
                iS_ri_hat_max_i - iS_ri_hat_min_i;
            R_0(i,i) = R_0(i,i) + iS_min_i + iS_max_i;

        end

        % k = 1,...,N-1
        row_min = u_cnstr.min_ineqRow_k;
        row_max = u_cnstr.max_ineqRow_k;

        for k = 1:N-1
            for i = 1:nu
                row_min_i = row_min(i);
                row_max_i = row_max(i);

                iS_min_i = iS_k(row_min_i,k);
                iS_max_i = iS_k(row_max_i,k);
                iS_ri_hat_min_i = iS_ri_hat_k(row_min_i,k);
                iS_ri_hat_max_i = iS_ri_hat_k(row_max_i,k);

                ru_k(i,k) = ru_k(i,k) + ...
                    iS_ri_hat_max_i - iS_ri_hat_min_i;
                R_k(i,i,k) = R_k(i,i,k) + iS_min_i + iS_max_i;

            end
        end
    elseif u_cnstr.min_limit
        
        % k = 0
        row = u_cnstr.min_ineqRow_0;

        for i = 1:nu
            row_i = row(i);

            iS_i = iS_0(row_i);
            iS_ri_hat_i = iS_ri_hat_0(row_i);

            % A_i = -I <-- u
            ru_0(i) = ru_0(i) - iS_ri_hat_i;
            R_0(i,i) =  R_0(i,i) + iS_i;

        end

        % k = 1,...,N-1
        row = u_cnstr.min_ineqRow_k;

        for k = 1:N-1
            for i = 1:nu

                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = -I <-- u
                ru_k(i,k) = ru_k(i,k) - iS_ri_hat_i;

                R_k(i,i,k) =  R_k(i,i,k) + iS_i;

            end
        end
    elseif u_cnstr.max_limit
        
        % k = 0
        row = u_cnstr.max_ineqRow_0;

        for i = 1:nu
            row_i = row(i);

            iS_i = iS_0(row_i);
            iS_ri_hat_i = iS_ri_hat_0(row_i);

            % A_i = I <-- u
            ru_0(i) = ru_0(i) + iS_ri_hat_i;
            R_0(i,i) =  R_0(i,i) + iS_i;

        end

        % k = 1,...,N-1
        row = u_cnstr.max_ineqRow_k;

        for k = 1:N-1
            for i = 1:nu

                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = I <-- u
                ru_k(i,k) = ru_k(i,k) + iS_ri_hat_i;

                R_k(i,i,k) =  R_k(i,i,k) + iS_i;

            end
        end
    end
end

if mpc.has_du_cnstr
    if du_cnstr.min_limit && du_cnstr.max_limit

        row_min = du_cnstr.min_ineqRow_0;
        row_max = du_cnstr.max_ineqRow_0;

        % k = 0
        for i = 1:nu
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            iS_min_i = iS_0(row_min_i);
            iS_max_i = iS_0(row_max_i);
            iS_ri_hat_min_i = iS_ri_hat_0(row_min_i);
            iS_ri_hat_max_i = iS_ri_hat_0(row_max_i);

            ru_0(i) = ru_0(i) + ...
                iS_ri_hat_max_i - iS_ri_hat_min_i;
            R_0(i,i) = R_0(i,i) + iS_min_i + iS_max_i;

        end

        % k = 1,...,N-1
        row_min = du_cnstr.min_ineqRow_k;
        row_max = du_cnstr.max_ineqRow_k;

        for k = 1:N-1
            for i = 1:nu
                su_i = su_col(i);
                row_min_i = row_min(i);
                row_max_i = row_max(i);

                iS_min_i = iS_k(row_min_i,k);
                iS_max_i = iS_k(row_max_i,k);
                iS_ri_hat_min_i = iS_ri_hat_k(row_min_i,k);
                iS_ri_hat_max_i = iS_ri_hat_k(row_max_i,k);

                iS_ri_hat_delta_i = iS_ri_hat_max_i - iS_ri_hat_min_i;
                iS_sum_i = iS_min_i + iS_max_i;

                rse_k(su_i,k) = rse_k(su_i,k) - ...
                    iS_ri_hat_delta_i;
                ru_k(i,k) = ru_k(i,k) + iS_ri_hat_delta_i;

                Q_k(su_i,su_i,k) = Q_k(su_i,su_i,k) + iS_sum_i;
                R_k(i,i,k) = R_k(i,i,k) + iS_sum_i;
                Y_k(i,su_i,k) = Y_k(i,su_i,k) - iS_sum_i;

            end
        end
    elseif du_cnstr.min_limit
        
        % k = 0
        row = du_cnstr.min_ineqRow_0;

        for i = 1:nu
            row_i = row(i);

            iS_i = iS_0(row_i);
            iS_ri_hat_i = iS_ri_hat_0(row_i);

            % A_i = -I <-- u
            ru_0(i) = ru_0(i) - iS_ri_hat_i;
            R_0(i,i) =  R_0(i,i) + iS_i;

        end

        % k = 1,...,N-1
        row = du_cnstr.min_ineqRow_k;

        for k = 1:N-1
            for i = 1:nu

                su_i = su_col(i);
                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = [I -I] <-- [su u]^T
                rse_k(su_i,k) = rse_k(su_i,k) + iS_ri_hat_i;
                ru_k(i,k) = ru_k(i,k) - iS_ri_hat_i;

                Q_k(su_i,su_i,k) = Q_k(su_i,su_i,k) + iS_i;
                R_k(i,i,k) =  R_k(i,i,k) + iS_i;
                Y_k(i,su_i,k) = Y_k(i,su_i,k) - iS_i;

            end
        end
    elseif du_cnstr.max_limit
        
        % k = 0
        row = du_cnstr.max_ineqRow_0;

        for i = 1:nu
            row_i = row(i);

            iS_i = iS_0(row_i);
            iS_ri_hat_i = iS_ri_hat_0(row_i);

            % A_i = I <-- u
            ru_0(i) = ru_0(i) + iS_ri_hat_i;
            R_0(i,i) =  R_0(i,i) + iS_i;

        end

        % k = 1,...,N-1
        row = du_cnstr.max_ineqRow_k;

        for k = 1:N-1
            for i = 1:nu

                su_i = su_col(i);
                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = [-I I] <-- [su u]^T
                rse_k(su_i,k) = rse_k(su_i,k) - iS_ri_hat_i;
                ru_k(i,k) = ru_k(i,k) + iS_ri_hat_i;

                Q_k(su_i,su_i,k) = Q_k(su_i,su_i,k) + iS_i;
                R_k(i,i,k) =  R_k(i,i,k) + iS_i;
                Y_k(i,su_i,k) = Y_k(i,su_i,k) - iS_i;

            end
        end
    end
end

if mpc.has_s_cnstr
    if s_cnstr.min_limit && s_cnstr.max_limit

        row_min = s_cnstr.min_ineqRow_k;
        row_max = s_cnstr.max_ineqRow_k;

        % k = 1,...,N-1
        for k = 1:N-1
            for i = 1:nx
                s_i = s_col(i);
                row_min_i = row_min(i);
                row_max_i = row_max(i);

                iS_min_i = iS_k(row_min_i,k);
                iS_max_i = iS_k(row_max_i,k);
                iS_ri_hat_min_i = iS_ri_hat_k(row_min_i,k);
                iS_ri_hat_max_i = iS_ri_hat_k(row_max_i,k);

                rse_k(s_i,k) = rse_k(s_i,k) + ...
                    iS_ri_hat_max_i - iS_ri_hat_min_i;
                Q_k(s_i,s_i,k) = Q_k(s_i,s_i,k) + ...
                    iS_min_i + iS_max_i;

            end
        end

        % k = N
        row_min = s_cnstr.min_ineqRow_ter;
        row_max = s_cnstr.max_ineqRow_ter;

        for i = 1:nx
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            iS_min_i = iS_ter(row_min_i);
            iS_max_i = iS_ter(row_max_i);
            iS_ri_hat_min_i = iS_ri_hat_ter(row_min_i);
            iS_ri_hat_max_i = iS_ri_hat_ter(row_max_i);

            rse_ter(i) = rse_ter(i) + ...
                iS_ri_hat_max_i - iS_ri_hat_min_i;
            Q_ter(i,i) = Q_ter(i,i) + iS_min_i + iS_max_i;

        end
    elseif s_cnstr.min_limit

        % k = 1,...,N-1
        row = s_cnstr.min_ineqRow_k;

        for k = 1:N-1
            for i = 1:nx

                s_i = s_col(i);
                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = -I <-- s
                rse_k(s_i,k) = rse_k(s_i,k) - iS_ri_hat_i;

                Q_k(s_i,s_i,k) = Q_k(s_i,s_i,k) + iS_i;

            end
        end

        % k = N
        row = s_cnstr.min_ineqRow_ter;

        for i = 1:nx
            row_i = row(i);

            iS_i = iS_ter(row_i);
            iS_ri_hat_i = iS_ri_hat_ter(row_i);

            % A_i = -I <-- s
            rse_ter(i) = rse_ter(i) - iS_ri_hat_i;
            Q_ter(i,i) =  Q_ter(i,i) + iS_i;

        end
    elseif s_cnstr.max_limit

        % k = 1,...,N-1
        row = s_cnstr.max_ineqRow_k;

        for k = 1:N-1
            for i = 1:nx

                s_i = s_col(i);
                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = I <-- s
                rse_k(s_i,k) = rse_k(s_i,k) + iS_ri_hat_i;

                Q_k(s_i,s_i,k) = Q_k(s_i,s_i,k) + iS_i;

            end
        end

        % k = N
        row = s_cnstr.max_ineqRow_ter;

        for i = 1:nx
            row_i = row(i);

            iS_i = iS_ter(row_i);
            iS_ri_hat_i = iS_ri_hat_ter(row_i);

            % A_i = I <-- s
            rse_ter(i) = rse_ter(i) + iS_ri_hat_i;
            Q_ter(i,i) =  Q_ter(i,i) + iS_i;

        end
    end
end


% if has_y_cnstr && y_use_k0
%     if y_cnstr.min_limit
%         row = y_cnstr.min_ineqRow_0;
%         ru_hat_0 = ru_hat_0 - D_0'*iS_ri_hat_0(row);
%     end
%     if y_cnstr.max_limit
%         row = y_cnstr.max_ineqRow_0;
%         ru_hat_0 = ru_hat_0 + D_0'*iS_ri_hat_0(row);
%     end
% end
% 
% if has_h_cnstr
%     if h_cnstr.min_limit && h_cnstr.use_k0
%         row = h_cnstr.min_ineqRow_0;
%         ru_hat_0 = ru_hat_0 - Dh_0'*iS_ri_hat_0(row);
%     end
%     if h_cnstr.max_limit && h_cnstr.use_k0
%         row = h_cnstr.max_ineqRow_0;
%         ru_hat_0 = ru_hat_0 + Dh_0'*iS_ri_hat_0(row);
%     end
% end



% 
% for k = 1:N-1
% 
%     rse_hat_k(:,k) = rse_k(:,k) + Ai_k(:,se_col,k)'*iS_ri_hat_k(:,k);
%     ru_hat_k(:,k) = ru_k(:,k) + Ai_k(:,u_col,k)'*iS_ri_hat_k(:,k);
% end

%rse_hat_ter = rse_ter + Ai_ter'*iS_ri_hat_ter;

%% H = Hess(f0) + Ai'*(S^-1)*Ai

% iS_Ai_0 = iS_0.*Ai_0;
% %R_0(:,:) = t*H_f0_0 + Ai_0'*iS_Ai_0;
% 
% for k = 1:N-1
% 
%     iS_Ai_k(:,:,k) = iS_k(:,k).*Ai_k(:,:,k);
% 
% %     H_k(:,:,k) = t*H_f0_k(:,:,k) + ...
% %                         Ai_k(:,:,k)'*iS_Ai_k(:,:,k);
% % 
% %     Q_k(:,:,k) = H_k(se_col,se_col,k);
% %     R_k(:,:,k) = H_k(u_col,u_col,k);
% %     Y_k(:,:,k) = H_k(u_col,se_col,k);
% end
% 
% %Q_ter(:,:) = t*H_f0_ter;
% if ng_k(3)
%     iS_Ai_ter = iS_ter.*Ai_ter;
%     %Q_ter(:,:) = Q_ter(:,:) + Ai_ter'*iS_Ai_ter;
% end

end
