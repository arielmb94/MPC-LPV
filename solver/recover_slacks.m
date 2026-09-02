function [delta_g_0,delta_g_k,delta_g_ter,...
          delta_v_0,delta_v_k,delta_v_ter] = recover_slacks(mpc,N,...
                    delta_g_0,delta_g_k,delta_g_ter,...
                    delta_v_0,delta_v_k,delta_v_ter,delta_se,delta_u,...
                    iS_ri_hat_0,iS_ri_hat_k,iS_ri_hat_ter,iS_0,iS_k,iS_ter,...
                    g_0,g_k,g_ter,g2_0,g2_k,g2_ter,v2_0,v2_k,v2_ter,...
                    rv_v2_0,rv_v2_k,rv_v2_ter,s_cnstr,u_cnstr,du_cnstr,...
                    nu,nx,su_col,y_cnstr,D_0,C,D,C_ter,s_col,...
                    h_cnstr,Dh_0,Ch,Dsuh,Dh,Ch_ter)
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

if mpc.has_y_cnstr
    if y_cnstr.min_limit && y_cnstr.max_limit

        % k = 0
        if y_cnstr.use_k0
            row_min = y_cnstr.min_ineqRow_0;
            row_max = y_cnstr.max_ineqRow_0;
            row_v_min = y_cnstr.min_row_v_0;
            row_v_max = y_cnstr.max_row_v_0;

            delta_y_0 = D_0 * delta_u(:,1);

            mu_min = iS_ri_hat_0(row_min) - iS_0(row_min) .* delta_y_0;
            mu_max = iS_ri_hat_0(row_max) + iS_0(row_max) .* delta_y_0;

            delta_g_0(row_min) = g_0(row_min) - g2_0(row_min) .* mu_min;
            delta_g_0(row_max) = g_0(row_max) - g2_0(row_max) .* mu_max;

            delta_v_0(row_v_min) = -rv_v2_0(row_v_min) + v2_0(row_v_min) .* mu_min;
            delta_v_0(row_v_max) = -rv_v2_0(row_v_max) + v2_0(row_v_max) .* mu_max;
        end

        % k = 1,...,N-1
        row_min = y_cnstr.min_ineqRow_k;
        row_max = y_cnstr.max_ineqRow_k;
        row_v_min = y_cnstr.min_row_v_k;
        row_v_max = y_cnstr.max_row_v_k;
        for k = 1:N-1
            if y_cnstr.use_s && y_cnstr.use_u
                delta_y_k = C(:,:,k) * delta_se(s_col,k) + D(:,:,k) * delta_u(:,k+1);
            elseif y_cnstr.use_s
                delta_y_k = C(:,:,k) * delta_se(s_col,k);
            elseif y_cnstr.use_u
                delta_y_k = D(:,:,k) * delta_u(:,k+1);
            else
                delta_y_k = 0;
            end

            mu_min = iS_ri_hat_k(row_min,k) - iS_k(row_min,k) .* delta_y_k;
            mu_max = iS_ri_hat_k(row_max,k) + iS_k(row_max,k) .* delta_y_k;

            delta_g_k(row_min,k) = g_k(row_min,k) - g2_k(row_min,k) .* mu_min;
            delta_g_k(row_max,k) = g_k(row_max,k) - g2_k(row_max,k) .* mu_max;

            delta_v_k(row_v_min,k) = -rv_v2_k(row_v_min,k) + v2_k(row_v_min,k) .* mu_min;
            delta_v_k(row_v_max,k) = -rv_v2_k(row_v_max,k) + v2_k(row_v_max,k) .* mu_max;
        end

        % k = N
        if y_cnstr.use_ter
            row_min = y_cnstr.min_ineqRow_ter;
            row_max = y_cnstr.max_ineqRow_ter;

            delta_y_ter = C_ter * delta_se(s_col,N);

            mu_min = iS_ri_hat_ter(row_min) - iS_ter(row_min) .* delta_y_ter;
            mu_max = iS_ri_hat_ter(row_max) + iS_ter(row_max) .* delta_y_ter;

            delta_g_ter(row_min) = g_ter(row_min) - g2_ter(row_min) .* mu_min;
            delta_g_ter(row_max) = g_ter(row_max) - g2_ter(row_max) .* mu_max;

            delta_v_ter(row_min) = -rv_v2_ter(row_min) + v2_ter(row_min) .* mu_min;
            delta_v_ter(row_max) = -rv_v2_ter(row_max) + v2_ter(row_max) .* mu_max;
        end

    elseif y_cnstr.min_limit

        % k = 0
        if y_cnstr.use_k0
            row = y_cnstr.min_ineqRow_0;
            row_v = y_cnstr.min_row_v_0;

            delta_y_0 = D_0 * delta_u(:,1);
            mu = iS_ri_hat_0(row) - iS_0(row) .* delta_y_0;

            delta_g_0(row) = g_0(row) - g2_0(row) .* mu;
            delta_v_0(row_v) = -rv_v2_0(row_v) + v2_0(row_v) .* mu;
        end

        % k = 1,...,N-1
        row = y_cnstr.min_ineqRow_k;
        row_v = y_cnstr.min_row_v_k;
        for k = 1:N-1
            if y_cnstr.use_s && y_cnstr.use_u
                delta_y_k = C(:,:,k) * delta_se(s_col,k) + D(:,:,k) * delta_u(:,k+1);
            elseif y_cnstr.use_s
                delta_y_k = C(:,:,k) * delta_se(s_col,k);
            elseif y_cnstr.use_u
                delta_y_k = D(:,:,k) * delta_u(:,k+1);
            else
                delta_y_k = 0;
            end

            mu = iS_ri_hat_k(row,k) - iS_k(row,k) .* delta_y_k;

            delta_g_k(row,k) = g_k(row,k) - g2_k(row,k) .* mu;
            delta_v_k(row_v,k) = -rv_v2_k(row_v,k) + v2_k(row_v,k) .* mu;
        end

        % k = N
        if y_cnstr.use_ter
            row = y_cnstr.min_ineqRow_ter;

            delta_y_ter = C_ter * delta_se(s_col,N);
            mu = iS_ri_hat_ter(row) - iS_ter(row) .* delta_y_ter;

            delta_g_ter(row) = g_ter(row) - g2_ter(row) .* mu;
            delta_v_ter(row) = -rv_v2_ter(row) + v2_ter(row) .* mu;
        end

    elseif y_cnstr.max_limit

        % k = 0
        if y_cnstr.use_k0
            row = y_cnstr.max_ineqRow_0;
            row_v = y_cnstr.max_row_v_0;

            delta_y_0 = D_0 * delta_u(:,1);
            mu = iS_ri_hat_0(row) + iS_0(row) .* delta_y_0;

            delta_g_0(row) = g_0(row) - g2_0(row) .* mu;
            delta_v_0(row_v) = -rv_v2_0(row_v) + v2_0(row_v) .* mu;
        end

        % k = 1,...,N-1
        row = y_cnstr.max_ineqRow_k;
        row_v = y_cnstr.max_row_v_k;
        for k = 1:N-1
            if y_cnstr.use_s && y_cnstr.use_u
                delta_y_k = C(:,:,k) * delta_se(s_col,k) + D(:,:,k) * delta_u(:,k+1);
            elseif y_cnstr.use_s
                delta_y_k = C(:,:,k) * delta_se(s_col,k);
            elseif y_cnstr.use_u
                delta_y_k = D(:,:,k) * delta_u(:,k+1);
            else
                delta_y_k = 0;
            end

            mu = iS_ri_hat_k(row,k) + iS_k(row,k) .* delta_y_k;

            delta_g_k(row,k) = g_k(row,k) - g2_k(row,k) .* mu;
            delta_v_k(row_v,k) = -rv_v2_k(row_v,k) + v2_k(row_v,k) .* mu;
        end

        % k = N
        if y_cnstr.use_ter
            row = y_cnstr.max_ineqRow_ter;

            delta_y_ter = C_ter * delta_se(s_col,N);
            mu = iS_ri_hat_ter(row) + iS_ter(row) .* delta_y_ter;

            delta_g_ter(row) = g_ter(row) - g2_ter(row) .* mu;
            delta_v_ter(row) = -rv_v2_ter(row) + v2_ter(row) .* mu;
        end

    end
end

if mpc.has_h_cnstr
    if h_cnstr.min_limit && h_cnstr.max_limit

        % k = 0
        if h_cnstr.use_k0
            row_min = h_cnstr.min_ineqRow_0;
            row_max = h_cnstr.max_ineqRow_0;
            row_v_min = h_cnstr.min_row_v_0;
            row_v_max = h_cnstr.max_row_v_0;

            delta_h_0 = Dh_0 * delta_u(:,1);

            mu_min = iS_ri_hat_0(row_min) - iS_0(row_min) .* delta_h_0;
            mu_max = iS_ri_hat_0(row_max) + iS_0(row_max) .* delta_h_0;

            delta_g_0(row_min) = g_0(row_min) - g2_0(row_min) .* mu_min;
            delta_g_0(row_max) = g_0(row_max) - g2_0(row_max) .* mu_max;

            delta_v_0(row_v_min) = -rv_v2_0(row_v_min) + v2_0(row_v_min) .* mu_min;
            delta_v_0(row_v_max) = -rv_v2_0(row_v_max) + v2_0(row_v_max) .* mu_max;
        end

        % k = 1,...,N-1
        row_min = h_cnstr.min_ineqRow_k;
        row_max = h_cnstr.max_ineqRow_k;
        row_v_min = h_cnstr.min_row_v_k;
        row_v_max = h_cnstr.max_row_v_k;
        for k = 1:N-1
            if h_cnstr.use_s && h_cnstr.use_su && h_cnstr.use_u
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + ...
                    Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
            elseif h_cnstr.use_s && h_cnstr.use_su
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dsuh(:,:,k) * delta_se(su_col,k);
            elseif h_cnstr.use_s && h_cnstr.use_u
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dh(:,:,k) * delta_u(:,k+1);
            elseif h_cnstr.use_su && h_cnstr.use_u
                delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
            elseif h_cnstr.use_s
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k);
            elseif h_cnstr.use_su
                delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k);
            elseif h_cnstr.use_u
                delta_h_k = Dh(:,:,k) * delta_u(:,k+1);
            end

            mu_min = iS_ri_hat_k(row_min,k) - iS_k(row_min,k) .* delta_h_k;
            mu_max = iS_ri_hat_k(row_max,k) + iS_k(row_max,k) .* delta_h_k;

            delta_g_k(row_min,k) = g_k(row_min,k) - g2_k(row_min,k) .* mu_min;
            delta_g_k(row_max,k) = g_k(row_max,k) - g2_k(row_max,k) .* mu_max;

            delta_v_k(row_v_min,k) = -rv_v2_k(row_v_min,k) + v2_k(row_v_min,k) .* mu_min;
            delta_v_k(row_v_max,k) = -rv_v2_k(row_v_max,k) + v2_k(row_v_max,k) .* mu_max;
        end

        % k = N
        if h_cnstr.use_ter
            row_min = h_cnstr.min_ineqRow_ter;
            row_max = h_cnstr.max_ineqRow_ter;

            delta_h_ter = Ch_ter * delta_se(s_col,N);

            mu_min = iS_ri_hat_ter(row_min) - iS_ter(row_min) .* delta_h_ter;
            mu_max = iS_ri_hat_ter(row_max) + iS_ter(row_max) .* delta_h_ter;

            delta_g_ter(row_min) = g_ter(row_min) - g2_ter(row_min) .* mu_min;
            delta_g_ter(row_max) = g_ter(row_max) - g2_ter(row_max) .* mu_max;

            delta_v_ter(row_min) = -rv_v2_ter(row_min) + v2_ter(row_min) .* mu_min;
            delta_v_ter(row_max) = -rv_v2_ter(row_max) + v2_ter(row_max) .* mu_max;
        end

    elseif h_cnstr.min_limit

        % k = 0
        if h_cnstr.use_k0
            row = h_cnstr.min_ineqRow_0;
            row_v = h_cnstr.min_row_v_0;

            delta_h_0 = Dh_0 * delta_u(:,1);
            mu = iS_ri_hat_0(row) - iS_0(row) .* delta_h_0;

            delta_g_0(row) = g_0(row) - g2_0(row) .* mu;
            delta_v_0(row_v) = -rv_v2_0(row_v) + v2_0(row_v) .* mu;
        end

        % k = 1,...,N-1
        row = h_cnstr.min_ineqRow_k;
        row_v = h_cnstr.min_row_v_k;
        for k = 1:N-1
            if h_cnstr.use_s && h_cnstr.use_su && h_cnstr.use_u
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + ...
                    Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
            elseif h_cnstr.use_s && h_cnstr.use_su
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dsuh(:,:,k) * delta_se(su_col,k);
            elseif h_cnstr.use_s && h_cnstr.use_u
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dh(:,:,k) * delta_u(:,k+1);
            elseif h_cnstr.use_su && h_cnstr.use_u
                delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
            elseif h_cnstr.use_s
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k);
            elseif h_cnstr.use_su
                delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k);
            elseif h_cnstr.use_u
                delta_h_k = Dh(:,:,k) * delta_u(:,k+1);
            end

            mu = iS_ri_hat_k(row,k) - iS_k(row,k) .* delta_h_k;

            delta_g_k(row,k) = g_k(row,k) - g2_k(row,k) .* mu;
            delta_v_k(row_v,k) = -rv_v2_k(row_v,k) + v2_k(row_v,k) .* mu;
        end

        % k = N
        if h_cnstr.use_ter
            row = h_cnstr.min_ineqRow_ter;

            delta_h_ter = Ch_ter * delta_se(s_col,N);
            mu = iS_ri_hat_ter(row) - iS_ter(row) .* delta_h_ter;

            delta_g_ter(row) = g_ter(row) - g2_ter(row) .* mu;
            delta_v_ter(row) = -rv_v2_ter(row) + v2_ter(row) .* mu;
        end

    elseif h_cnstr.max_limit

        % k = 0
        if h_cnstr.use_k0
            row = h_cnstr.max_ineqRow_0;
            row_v = h_cnstr.max_row_v_0;

            delta_h_0 = Dh_0 * delta_u(:,1);
            mu = iS_ri_hat_0(row) + iS_0(row) .* delta_h_0;

            delta_g_0(row) = g_0(row) - g2_0(row) .* mu;
            delta_v_0(row_v) = -rv_v2_0(row_v) + v2_0(row_v) .* mu;
        end

        % k = 1,...,N-1
        row = h_cnstr.max_ineqRow_k;
        row_v = h_cnstr.max_row_v_k;
        for k = 1:N-1
            if h_cnstr.use_s && h_cnstr.use_su && h_cnstr.use_u
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + ...
                    Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
            elseif h_cnstr.use_s && h_cnstr.use_su
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dsuh(:,:,k) * delta_se(su_col,k);
            elseif h_cnstr.use_s && h_cnstr.use_u
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k) + Dh(:,:,k) * delta_u(:,k+1);
            elseif h_cnstr.use_su && h_cnstr.use_u
                delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k) + Dh(:,:,k) * delta_u(:,k+1);
            elseif h_cnstr.use_s
                delta_h_k = Ch(:,:,k) * delta_se(s_col,k);
            elseif h_cnstr.use_su
                delta_h_k = Dsuh(:,:,k) * delta_se(su_col,k);
            elseif h_cnstr.use_u
                delta_h_k = Dh(:,:,k) * delta_u(:,k+1);
            end

            mu = iS_ri_hat_k(row,k) + iS_k(row,k) .* delta_h_k;

            delta_g_k(row,k) = g_k(row,k) - g2_k(row,k) .* mu;
            delta_v_k(row_v,k) = -rv_v2_k(row_v,k) + v2_k(row_v,k) .* mu;
        end

        % k = N
        if h_cnstr.use_ter
            row = h_cnstr.max_ineqRow_ter;

            delta_h_ter = Ch_ter * delta_se(s_col,N);
            mu = iS_ri_hat_ter(row) + iS_ter(row) .* delta_h_ter;

            delta_g_ter(row) = g_ter(row) - g2_ter(row) .* mu;
            delta_v_ter(row) = -rv_v2_ter(row) + v2_ter(row) .* mu;
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
