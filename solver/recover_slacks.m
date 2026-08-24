function mpc = recover_slacks(mpc,delta_u,delta_se)
% Cache read-only arrays and output workspaces for the hot loops.
delta_g_0 = mpc.delta_g_0;
delta_g_k = mpc.delta_g_k;
delta_g_ter = mpc.delta_g_ter;
delta_v_0 = mpc.delta_v_0;
delta_v_k = mpc.delta_v_k;
delta_v_ter = mpc.delta_v_ter;
if mpc.ng_k(1)
    iS_ri_hat_0 = mpc.iS_ri_hat_0;
    g_0 = mpc.g_0;
    g2_0 = mpc.g2_0;
end
if mpc.ng_k(2)
    iS_k = mpc.iS_k;
    iS_ri_hat_k = mpc.iS_ri_hat_k;
    g_k = mpc.g_k;
    g2_k = mpc.g2_k;
end
if mpc.ng_k(3)
    iS_ter = mpc.iS_ter;
    iS_ri_hat_ter = mpc.iS_ri_hat_ter;
    g_ter = mpc.g_ter;
    g2_ter = mpc.g2_ter;
end
v2_k = mpc.v2_k;
v2_ter = mpc.v2_ter;
rv_v2_k = mpc.rv_v2_k;
rv_v2_ter = mpc.rv_v2_ter;
su_col = mpc.su_col;

%% \mu_i = (-S)^{-1}(-\hat{r}_i-A_i\Delta x) = S^{-1}(\hat{r}_i+A_i\Delta x)
% Delta g  = (g^2)(-(-1/g)-\mu_i) = g - g^2*mu
% \Delta v = (v^2)(-r_v+ \mu_i) = -v^2*r_v + v^2*mu_i

%k = 0
%mpc.mu_i_0(:) = mpc.iS_ri_hat_0+mpc.iS_Ai_0*delta_u(:,1);

if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit && mpc.s_cnstr.max_limit

        %k = 1,...N-1
        row_min = mpc.s_cnstr.min_ineqRow_k;
        row_max = mpc.s_cnstr.max_ineqRow_k;
        row_v_min = mpc.s_cnstr.min_row_v_k;
        row_v_max = mpc.s_cnstr.max_row_v_k;
        for k = 1:mpc.N-1
            for i = 1:mpc.nx
                row_min_i = row_min(i);
                row_max_i = row_max(i);

                mu_min_i = iS_ri_hat_k(row_min_i,k) - ...
                    iS_k(row_min_i,k)*delta_se(i,k);
                mu_max_i = iS_ri_hat_k(row_max_i,k) + ...
                    iS_k(row_max_i,k)*delta_se(i,k);

                g_min_i = g_k(row_min_i,k);
                g_max_i = g_k(row_max_i,k);
                g2_min_i = g2_k(row_min_i,k);
                g2_max_i = g2_k(row_max_i,k);

                delta_g_k(row_min_i,k) = g_min_i - g2_min_i*mu_min_i;
                delta_g_k(row_max_i,k) = g_max_i - g2_max_i*mu_max_i;

                row_v_min_i = row_v_min(i);
                row_v_max_i = row_v_max(i);
                v2_min_i = v2_k(row_v_min_i,k);
                v2_max_i = v2_k(row_v_max_i,k);

                delta_v_k(row_v_min_i,k) = ...
                    -rv_v2_k(row_v_min_i,k) + v2_min_i*mu_min_i;
                delta_v_k(row_v_max_i,k) = ...
                    -rv_v2_k(row_v_max_i,k) + v2_max_i*mu_max_i;
            end
        end

        %k = N
        row_min = mpc.s_cnstr.min_ineqRow_ter;
        row_max = mpc.s_cnstr.max_ineqRow_ter;
        for i = 1:mpc.nx
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            mu_min_i = iS_ri_hat_ter(row_min_i) - ...
                iS_ter(row_min_i)*delta_se(i,mpc.N);
            mu_max_i = iS_ri_hat_ter(row_max_i) + ...
                iS_ter(row_max_i)*delta_se(i,mpc.N);

            g_min_i = g_ter(row_min_i);
            g_max_i = g_ter(row_max_i);
            g2_min_i = g2_ter(row_min_i);
            g2_max_i = g2_ter(row_max_i);

            delta_g_ter(row_min_i) = g_min_i - g2_min_i*mu_min_i;
            delta_g_ter(row_max_i) = g_max_i - g2_max_i*mu_max_i;

            v2_min_i = v2_ter(row_min_i);
            v2_max_i = v2_ter(row_max_i);

            delta_v_ter(row_min_i) = ...
                -rv_v2_ter(row_min_i) + v2_min_i*mu_min_i;
            delta_v_ter(row_max_i) = ...
                -rv_v2_ter(row_max_i) + v2_max_i*mu_max_i;
        end
    elseif mpc.s_cnstr.min_limit
        
        %k = 1,...N-1
        row = mpc.s_cnstr.min_ineqRow_k;
        row_v = mpc.s_cnstr.min_row_v_k;
        for k = 1:mpc.N-1
            for i = 1:mpc.nx
                row_i = row(i);
                
                % A_i = -I <-- s
                mu_i = iS_ri_hat_k(row_i,k)-iS_k(row_i,k)*delta_se(i,k);

                g_i = g_k(row_i,k);
                g2_i = g2_k(row_i,k);

                delta_g_k(row_i,k) = g_i - g2_i*mu_i;

                row_v_i = row_v(i);
                v2_i = v2_k(row_v_i,k);
                
                delta_v_k(row_v_i,k) = -rv_v2_k(row_v_i,k) + v2_i*mu_i;
            end
        end

        %k = N
        row = mpc.s_cnstr.min_ineqRow_ter;
        for i = 1:mpc.nx
            row_i = row(i);
            
            % A_i = -I <-- s
            mu_i = iS_ri_hat_ter(row_i)-iS_ter(row_i)*delta_se(i,mpc.N);

            g_i = g_ter(row_i);
            g2_i = g2_ter(row_i);
            delta_g_ter(row_i) = g_i - g2_i*mu_i;

            v2_i = v2_ter(row_i);
            delta_v_ter(row_i) = -rv_v2_ter(row_i) + v2_i*mu_i;
        end
    elseif mpc.s_cnstr.max_limit
        
        %k = 1,...N-1
        row = mpc.s_cnstr.max_ineqRow_k;
        row_v = mpc.s_cnstr.max_row_v_k;
        for k = 1:mpc.N-1
            for i = 1:mpc.nx
                row_i = row(i);

                % A_i = I <-- s
                mu_i = iS_ri_hat_k(row_i,k)+iS_k(row_i,k)*delta_se(i,k);

                g_i = g_k(row_i,k);
                g2_i = g2_k(row_i,k);
                delta_g_k(row_i,k) = g_i - g2_i*mu_i;

                row_v_i = row_v(i);
                v2_i = v2_k(row_v_i,k);

                delta_v_k(row_v_i,k) = -rv_v2_k(row_v_i,k) + v2_i*mu_i;
            end
        end

        %k = N
        row = mpc.s_cnstr.max_ineqRow_ter;
        for i = 1:mpc.nx
            row_i = row(i);

            % A_i = I <-- s
            mu_i = iS_ri_hat_ter(row_i)+iS_ter(row_i)*delta_se(i,mpc.N);

            g_i = g_ter(row_i);
            g2_i = g2_ter(row_i);
            delta_g_ter(row_i) = g_i - g2_i*mu_i;

            v2_i = v2_ter(row_i);
            delta_v_ter(row_i) = -rv_v2_ter(row_i) + v2_i*mu_i;
        end
    end
end

if mpc.has_u_cnstr
    if mpc.u_cnstr.min_limit && mpc.u_cnstr.max_limit
        %k = 0
        row_min = mpc.u_cnstr.min_ineqRow_0;
        row_max = mpc.u_cnstr.max_ineqRow_0;
        for i = 1:mpc.nu
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            g_min_i = g_0(row_min_i);
            g_max_i = g_0(row_max_i);
            g2_min_i = g2_0(row_min_i);
            g2_max_i = g2_0(row_max_i);

            delta_g_0(row_min_i) = g_min_i - ...
                g2_min_i*iS_ri_hat_0(row_min_i) + delta_u(i,1);
            delta_g_0(row_max_i) = g_max_i - ...
                g2_max_i*iS_ri_hat_0(row_max_i) - delta_u(i,1);
        end

        %k = 1,...N-1
        row_min = mpc.u_cnstr.min_ineqRow_k;
        row_max = mpc.u_cnstr.max_ineqRow_k;
        for k = 1:mpc.N-1
            for i = 1:mpc.nu
                row_min_i = row_min(i);
                row_max_i = row_max(i);

                g_min_i = g_k(row_min_i,k);
                g_max_i = g_k(row_max_i,k);
                g2_min_i = g2_k(row_min_i,k);
                g2_max_i = g2_k(row_max_i,k);

                delta_g_k(row_min_i,k) = g_min_i - ...
                    g2_min_i*iS_ri_hat_k(row_min_i,k) + delta_u(i,k+1);
                delta_g_k(row_max_i,k) = g_max_i - ...
                    g2_max_i*iS_ri_hat_k(row_max_i,k) - delta_u(i,k+1);
            end
        end
    elseif mpc.u_cnstr.min_limit
        %k = 0
        row = mpc.u_cnstr.min_ineqRow_0;
        for i = 1:mpc.nu
            row_i = row(i);

            % A_i = -I <-- u
            g_i = g_0(row_i);
            g2_i = g2_0(row_i);
            delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) + delta_u(i,1);
        end

        %k = 1,...N-1
        row = mpc.u_cnstr.min_ineqRow_k;
        for k = 1:mpc.N-1
            for i = 1:mpc.nu
                row_i = row(i);

                % A_i = -I <-- u
                g_i = g_k(row_i,k);
                g2_i = g2_k(row_i,k);
                delta_g_k(row_i,k) = g_i - ...
                    g2_i*iS_ri_hat_k(row_i,k) + delta_u(i,k+1);
            end
        end
    elseif mpc.u_cnstr.max_limit
        %k = 0
        row = mpc.u_cnstr.max_ineqRow_0;
        for i = 1:mpc.nu
            row_i = row(i);

            % A_i = I <-- u
            g_i = g_0(row_i);
            g2_i = g2_0(row_i);
            delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) - delta_u(i,1);
        end
        %k = 1,...N-1
        row = mpc.u_cnstr.max_ineqRow_k;
        for k = 1:mpc.N-1
            for i = 1:mpc.nu
                row_i = row(i);

                % A_i = I <-- u
                g_i = g_k(row_i,k);
                g2_i = g2_k(row_i,k);
                delta_g_k(row_i,k) = g_i - ...
                    g2_i*iS_ri_hat_k(row_i,k) - delta_u(i,k+1);
            end
        end
    end
end

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit && mpc.du_cnstr.max_limit
        %k = 0
        row_min = mpc.du_cnstr.min_ineqRow_0;
        row_max = mpc.du_cnstr.max_ineqRow_0;
        for i = 1:mpc.nu
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            g_min_i = g_0(row_min_i);
            g_max_i = g_0(row_max_i);
            g2_min_i = g2_0(row_min_i);
            g2_max_i = g2_0(row_max_i);

            delta_g_0(row_min_i) = g_min_i - ...
                g2_min_i*iS_ri_hat_0(row_min_i) + delta_u(i,1);
            delta_g_0(row_max_i) = g_max_i - ...
                g2_max_i*iS_ri_hat_0(row_max_i) - delta_u(i,1);
        end

        %k = 1,...N-1
        row_min = mpc.du_cnstr.min_ineqRow_k;
        row_max = mpc.du_cnstr.max_ineqRow_k;
        for k = 1:mpc.N-1
            for i = 1:mpc.nu
                row_min_i = row_min(i);
                row_max_i = row_max(i);

                su_i = su_col(i);
                du_i = delta_u(i,k+1)-delta_se(su_i,k);

                g_min_i = g_k(row_min_i,k);
                g_max_i = g_k(row_max_i,k);
                g2_min_i = g2_k(row_min_i,k);
                g2_max_i = g2_k(row_max_i,k);

                delta_g_k(row_min_i,k) = g_min_i - ...
                    g2_min_i*iS_ri_hat_k(row_min_i,k) + du_i;
                delta_g_k(row_max_i,k) = g_max_i - ...
                    g2_max_i*iS_ri_hat_k(row_max_i,k) - du_i;
            end
        end
    elseif mpc.du_cnstr.min_limit
        %k = 0
        row = mpc.du_cnstr.min_ineqRow_0;
        for i = 1:mpc.nu
            row_i = row(i);

            % A_i = -I <-- u
            g_i = g_0(row_i);
            g2_i = g2_0(row_i);
            delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) + delta_u(i,1);
        end

        %k = 1,...N-1
        row = mpc.du_cnstr.min_ineqRow_k;
        for k = 1:mpc.N-1
            for i = 1:mpc.nu
                row_i = row(i);

                % A_i = [I -I] <-- [su u]
                su_i = su_col(i);
                du_i = delta_u(i,k+1)-delta_se(su_i,k);

                g_i = g_k(row_i,k);
                g2_i = g2_k(row_i,k);
                delta_g_k(row_i,k) = g_i - ...
                    g2_i*iS_ri_hat_k(row_i,k) + du_i;
            end
        end
    elseif mpc.du_cnstr.max_limit
        %k = 0
        row = mpc.du_cnstr.max_ineqRow_0;
        for i = 1:mpc.nu
            row_i = row(i);

            % A_i = I <-- u
            g_i = g_0(row_i);
            g2_i = g2_0(row_i);
            delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) - delta_u(i,1);
        end
        %k = 1,...N-1
        row = mpc.du_cnstr.max_ineqRow_k;
        for k = 1:mpc.N-1
            for i = 1:mpc.nu
                row_i = row(i);

                % A_i = [-I I] <-- [su u]
                su_i = su_col(i);
                du_i = delta_u(i,k+1)-delta_se(su_i,k);

                g_i = g_k(row_i,k);
                g2_i = g2_k(row_i,k);
                delta_g_k(row_i,k) = g_i - ...
                    g2_i*iS_ri_hat_k(row_i,k) - du_i;
            end
        end
    end
end


mpc.delta_g_0 = delta_g_0;
mpc.delta_g_k = delta_g_k;
mpc.delta_g_ter = delta_g_ter;
mpc.delta_v_0 = delta_v_0;
mpc.delta_v_k = delta_v_k;
mpc.delta_v_ter = delta_v_ter;

% for k = 1:mpc.N-1
%     mpc.mu_i_k(:,k) = mpc.iS_ri_hat_k(:,k) +...
%                       mpc.iS_Ai_k(:,mpc.se_col,k)*delta_se(:,k) +...
%                       mpc.iS_Ai_k(:,mpc.u_col,k)*delta_u(:,k+1);
% end
% %k = N
% if mpc.ng_k(3)
%     mpc.mu_i_ter(:) = mpc.iS_ri_hat_ter+mpc.iS_Ai_ter*delta_se(:,mpc.N);
% end
% 
% %% \Delta g  = (g^2)(-r_g-\mu_i) 
% %  \Delta v = (v^2)(-r_v+ \mu_i) 
% 
% mpc.delta_g_0(:) = mpc.g_0 - (mpc.g_0.^2).*mpc.mu_i_0;
% if mpc.nv_k(1)
%     mpc.delta_v_0(:) = -mpc.rv_v2_0+(mpc.v_0.^2).*mpc.mu_i_0(mpc.v_rows_0);
% end
% 
% mpc.delta_g_k = mpc.g_k - (mpc.g_k.^2).*mpc.mu_i_k;
% if mpc.nv_k(2)
%     mpc.delta_v_k(:) = -mpc.rv_v2_k+(mpc.v_k.^2).*mpc.mu_i_k(mpc.v_rows_k,:);
% end
% 
% if mpc.ng_k(3)
%     mpc.delta_g_ter = mpc.g_ter - (mpc.g_ter.^2).*mpc.mu_i_ter;
%     mpc.delta_v_ter = -mpc.rv_v2_ter +(mpc.v_ter.^2).*mpc.mu_i_ter;
% end

end
