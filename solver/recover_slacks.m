function [delta_g_0,delta_g_k,delta_g_ter,...
          delta_v_0,delta_v_k,delta_v_ter] = recover_slacks(mpc,N,...
                    delta_g_0,delta_g_k,delta_g_ter,...
                    delta_v_0,delta_v_k,delta_v_ter,delta_se,delta_u,...
                    iS_ri_hat_0,iS_ri_hat_k,iS_ri_hat_ter,iS_0,iS_k,iS_ter,...
                    g_0,g_k,g_ter,g2_0,g2_k,g2_ter,v2_0,v2_k,v2_ter,...
                    rv_v2_0,rv_v2_k,rv_v2_ter,s_cnstr,u_cnstr,du_cnstr,...
                    nu,nx,su_col)
%% \mu_i = (-S)^{-1}(-\hat{r}_i-A_i\Delta x) = S^{-1}(\hat{r}_i+A_i\Delta x)
% Delta g  = (g^2)(-(-1/g)-\mu_i) = g - g^2*mu
% \Delta v = (v^2)(-r_v+ \mu_i) = -v^2*r_v + v^2*mu_i

%k = 0
%mu_i_0(:) = iS_ri_hat_0+iS_Ai_0*delta_u(:,1);

if mpc.has_s_cnstr
    if s_cnstr.min_limit && s_cnstr.max_limit

        %k = 1,...N-1
        row_min = s_cnstr.min_ineqRow_k;
        row_max = s_cnstr.max_ineqRow_k;
        row_v_min = s_cnstr.min_row_v_k;
        row_v_max = s_cnstr.max_row_v_k;
        for k = 1:N-1
            for i = 1:nx
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
        row_min = s_cnstr.min_ineqRow_ter;
        row_max = s_cnstr.max_ineqRow_ter;
        for i = 1:nx
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            mu_min_i = iS_ri_hat_ter(row_min_i) - ...
                iS_ter(row_min_i)*delta_se(i,N);
            mu_max_i = iS_ri_hat_ter(row_max_i) + ...
                iS_ter(row_max_i)*delta_se(i,N);

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
    elseif s_cnstr.min_limit
        
        %k = 1,...N-1
        row = s_cnstr.min_ineqRow_k;
        row_v = s_cnstr.min_row_v_k;
        for k = 1:N-1
            for i = 1:nx
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
        row = s_cnstr.min_ineqRow_ter;
        for i = 1:nx
            row_i = row(i);
            
            % A_i = -I <-- s
            mu_i = iS_ri_hat_ter(row_i)-iS_ter(row_i)*delta_se(i,N);

            g_i = g_ter(row_i);
            g2_i = g2_ter(row_i);
            delta_g_ter(row_i) = g_i - g2_i*mu_i;

            v2_i = v2_ter(row_i);
            delta_v_ter(row_i) = -rv_v2_ter(row_i) + v2_i*mu_i;
        end
    elseif s_cnstr.max_limit
        
        %k = 1,...N-1
        row = s_cnstr.max_ineqRow_k;
        row_v = s_cnstr.max_row_v_k;
        for k = 1:N-1
            for i = 1:nx
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
        row = s_cnstr.max_ineqRow_ter;
        for i = 1:nx
            row_i = row(i);

            % A_i = I <-- s
            mu_i = iS_ri_hat_ter(row_i)+iS_ter(row_i)*delta_se(i,N);

            g_i = g_ter(row_i);
            g2_i = g2_ter(row_i);
            delta_g_ter(row_i) = g_i - g2_i*mu_i;

            v2_i = v2_ter(row_i);
            delta_v_ter(row_i) = -rv_v2_ter(row_i) + v2_i*mu_i;
        end
    end
end

if mpc.has_u_cnstr
    if u_cnstr.min_limit && u_cnstr.max_limit
        %k = 0
        row_min = u_cnstr.min_ineqRow_0;
        row_max = u_cnstr.max_ineqRow_0;
        for i = 1:nu
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
        row_min = u_cnstr.min_ineqRow_k;
        row_max = u_cnstr.max_ineqRow_k;
        for k = 1:N-1
            for i = 1:nu
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
    elseif u_cnstr.min_limit
        %k = 0
        row = u_cnstr.min_ineqRow_0;
        for i = 1:nu
            row_i = row(i);

            % A_i = -I <-- u
            g_i = g_0(row_i);
            g2_i = g2_0(row_i);
            delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) + delta_u(i,1);
        end

        %k = 1,...N-1
        row = u_cnstr.min_ineqRow_k;
        for k = 1:N-1
            for i = 1:nu
                row_i = row(i);

                % A_i = -I <-- u
                g_i = g_k(row_i,k);
                g2_i = g2_k(row_i,k);
                delta_g_k(row_i,k) = g_i - ...
                    g2_i*iS_ri_hat_k(row_i,k) + delta_u(i,k+1);
            end
        end
    elseif u_cnstr.max_limit
        %k = 0
        row = u_cnstr.max_ineqRow_0;
        for i = 1:nu
            row_i = row(i);

            % A_i = I <-- u
            g_i = g_0(row_i);
            g2_i = g2_0(row_i);
            delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) - delta_u(i,1);
        end
        %k = 1,...N-1
        row = u_cnstr.max_ineqRow_k;
        for k = 1:N-1
            for i = 1:nu
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
    if du_cnstr.min_limit && du_cnstr.max_limit
        %k = 0
        row_min = du_cnstr.min_ineqRow_0;
        row_max = du_cnstr.max_ineqRow_0;
        for i = 1:nu
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
        row_min = du_cnstr.min_ineqRow_k;
        row_max = du_cnstr.max_ineqRow_k;
        for k = 1:N-1
            for i = 1:nu
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
    elseif du_cnstr.min_limit
        %k = 0
        row = du_cnstr.min_ineqRow_0;
        for i = 1:nu
            row_i = row(i);

            % A_i = -I <-- u
            g_i = g_0(row_i);
            g2_i = g2_0(row_i);
            delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) + delta_u(i,1);
        end

        %k = 1,...N-1
        row = du_cnstr.min_ineqRow_k;
        for k = 1:N-1
            for i = 1:nu
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
    elseif du_cnstr.max_limit
        %k = 0
        row = du_cnstr.max_ineqRow_0;
        for i = 1:nu
            row_i = row(i);

            % A_i = I <-- u
            g_i = g_0(row_i);
            g2_i = g2_0(row_i);
            delta_g_0(row_i) = g_i - g2_i*iS_ri_hat_0(row_i) - delta_u(i,1);
        end
        %k = 1,...N-1
        row = du_cnstr.max_ineqRow_k;
        for k = 1:N-1
            for i = 1:nu
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


% for k = 1:N-1
%     mu_i_k(:,k) = iS_ri_hat_k(:,k) +...
%                       iS_Ai_k(:,se_col,k)*delta_se(:,k) +...
%                       iS_Ai_k(:,u_col,k)*delta_u(:,k+1);
% end
% %k = N
% if ng_k(3)
%     mu_i_ter(:) = iS_ri_hat_ter+iS_Ai_ter*delta_se(:,N);
% end
% 
% %% \Delta g  = (g^2)(-r_g-\mu_i) 
% %  \Delta v = (v^2)(-r_v+ \mu_i) 
% 
% delta_g_0(:) = g_0 - (g_0.^2).*mu_i_0;
% if nv_k(1)
%     delta_v_0(:) = -rv_v2_0+(v_0.^2).*mu_i_0(v_rows_0);
% end
% 
% delta_g_k = g_k - (g_k.^2).*mu_i_k;
% if nv_k(2)
%     delta_v_k(:) = -rv_v2_k+(v_k.^2).*mu_i_k(v_rows_k,:);
% end
% 
% if ng_k(3)
%     delta_g_ter = g_ter - (g_ter.^2).*mu_i_ter;
%     delta_v_ter = -rv_v2_ter +(v_ter.^2).*mu_i_ter;
% end

end
