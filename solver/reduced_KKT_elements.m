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
                        iS_0,iS_k,iS_ter,iS_ri_hat_0,iS_ri_hat_k,iS_ri_hat_ter,...
                        y_cnstr,D_0,C,D,C_ter,...
                        h_cnstr,Dh_0,Ch,Dsuh,Dh,Ch_ter)

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

%% rx_hat = t*grad_x + Ai'*S^-1*ri_hat

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

if mpc.has_y_cnstr

    if y_cnstr.min_limit && y_cnstr.max_limit
        if y_cnstr.use_k0
            row_min = y_cnstr.min_ineqRow_0;
            row_max = y_cnstr.max_ineqRow_0;

            iS_ri_hat_y0 = iS_ri_hat_0(row_max) - iS_ri_hat_0(row_min);

            iS_y0 = iS_0(row_min) + iS_0(row_max);

            ru_0 = ru_0 + D_0' * iS_ri_hat_y0;
            R_0  = R_0  + D_0' * (iS_y0 .* D_0);
        end

        % k = 1,...,N-1
        row_min = y_cnstr.min_ineqRow_k;
        row_max = y_cnstr.max_ineqRow_k;
        for k = 1:N-1

            iS_ri_hat_yk = iS_ri_hat_k(row_max,k)-iS_ri_hat_k(row_min,k);

            iS_yk = iS_k(row_min,k) + iS_k(row_max,k);

            if y_cnstr.use_s && y_cnstr.use_u

                rse_k(s_col,k) = rse_k(s_col,k) + C(:,:,k)' * iS_ri_hat_yk;
                ru_k(:,k) = ru_k(:,k) + D(:,:,k)' * iS_ri_hat_yk;

                iS_C_k = iS_yk .* C(:,:,k);

                Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k)  + C(:,:,k)' * iS_C_k;
                R_k(:,:,k)  = R_k(:,:,k)  + D(:,:,k)' * (iS_yk .* D(:,:,k));
                Y_k(:,s_col,k)  = Y_k(:,s_col,k)  + D(:,:,k)' * iS_C_k;
            elseif y_cnstr.use_s

                rse_k(s_col,k) = rse_k(s_col,k) + C(:,:,k)' * iS_ri_hat_yk;

                Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k)  + C(:,:,k)' * (iS_yk .* C(:,:,k));
            elseif y_cnstr.use_u

                ru_k(:,k) = ru_k(:,k) + D(:,:,k)' * iS_ri_hat_yk;

                R_k(:,:,k)  = R_k(:,:,k)  + D(:,:,k)' * (iS_yk .* D(:,:,k));
            end
        end

        % k = N
        if y_cnstr.use_ter
            row_min = y_cnstr.min_ineqRow_ter;
            row_max = y_cnstr.max_ineqRow_ter;

            iS_ri_hat_yter = iS_ri_hat_ter(row_max) - iS_ri_hat_ter(row_min);

            iS_yter = iS_ter(row_min) + iS_ter(row_max);

            rse_ter(s_col) = rse_ter(s_col) + C_ter' * iS_ri_hat_yter;
            Q_ter(s_col,s_col)  = Q_ter(s_col,s_col) + C_ter' * (iS_yter .* C_ter);
        end

    elseif y_cnstr.min_limit

        if y_cnstr.use_k0
            row_min = y_cnstr.min_ineqRow_0;

            ru_0 = ru_0 - D_0' * iS_ri_hat_0(row_min);
            R_0  = R_0  + D_0' * (iS_0(row_min) .* D_0);
        end

        % k = 1,...,N-1
        row_min = y_cnstr.min_ineqRow_k;
        for k = 1:N-1

            if y_cnstr.use_s && y_cnstr.use_u

                rse_k(s_col,k) = rse_k(s_col,k) - C(:,:,k)' * iS_ri_hat_k(row_min,k);
                ru_k(:,k) = ru_k(:,k) - D(:,:,k)' * iS_ri_hat_k(row_min,k);

                iS_C_k = iS_k(row_min,k) .* C(:,:,k);

                Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k) + C(:,:,k)' * iS_C_k;
                R_k(:,:,k)  = R_k(:,:,k) + D(:,:,k)' * (iS_k(row_min,k) .* D(:,:,k));
                Y_k(:,s_col,k)  = Y_k(:,s_col,k) + D(:,:,k)' * iS_C_k;
            elseif y_cnstr.use_s

                rse_k(s_col,k) = rse_k(s_col,k) - C(:,:,k)' * iS_ri_hat_k(row_min,k);

                Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k) + C(:,:,k)' * (iS_k(row_min,k) .* C(:,:,k));
            elseif y_cnstr.use_u

                ru_k(:,k) = ru_k(:,k) - D(:,:,k)' * iS_ri_hat_k(row_min,k);

                R_k(:,:,k) = R_k(:,:,k) + D(:,:,k)' * (iS_k(row_min,k) .* D(:,:,k));
            end
        end

        % k = N
        if y_cnstr.use_ter
            row_min = y_cnstr.min_ineqRow_ter;

            rse_ter(s_col) = rse_ter(s_col) - C_ter' * iS_ri_hat_ter(row_min);
            Q_ter(s_col,s_col)  = Q_ter(s_col,s_col) + C_ter' * (iS_ter(row_min) .* C_ter);
        end

    elseif y_cnstr.max_limit

        if y_cnstr.use_k0
            row_max = y_cnstr.max_ineqRow_0;

            ru_0 = ru_0 + D_0' * iS_ri_hat_0(row_max);
            R_0  = R_0  + D_0' * (iS_0(row_max) .* D_0);
        end

        % k = 1,...,N-1
        row_max = y_cnstr.max_ineqRow_k;
        for k = 1:N-1

            if y_cnstr.use_s && y_cnstr.use_u

                rse_k(s_col,k) = rse_k(s_col,k) + C(:,:,k)' * iS_ri_hat_k(row_max,k);
                ru_k(:,k) = ru_k(:,k) + D(:,:,k)' * iS_ri_hat_k(row_max,k);

                iS_C_k = iS_k(row_max,k) .* C(:,:,k);

                Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k) + C(:,:,k)' * iS_C_k;
                R_k(:,:,k)  = R_k(:,:,k) + D(:,:,k)' * (iS_k(row_max,k) .* D(:,:,k));
                Y_k(:,s_col,k)  = Y_k(:,s_col,k) + D(:,:,k)' * iS_C_k;
            elseif y_cnstr.use_s

                rse_k(s_col,k) = rse_k(s_col,k) + C(:,:,k)' * iS_ri_hat_k(row_max,k);

                Q_k(s_col,s_col,k)  = Q_k(s_col,s_col,k) + C(:,:,k)' * (iS_k(row_max,k) .* C(:,:,k));
            elseif y_cnstr.use_u

                ru_k(:,k) = ru_k(:,k) + D(:,:,k)' * iS_ri_hat_k(row_max,k);

                R_k(:,:,k) = R_k(:,:,k) + D(:,:,k)' * (iS_k(row_max,k) .* D(:,:,k));
            end
        end

        % k = N
        if y_cnstr.use_ter
            row_max = y_cnstr.max_ineqRow_ter;

            rse_ter(s_col) = rse_ter(s_col) + C_ter' * iS_ri_hat_ter(row_max);
            Q_ter(s_col,s_col)  = Q_ter(s_col,s_col) + C_ter' * (iS_ter(row_max) .* C_ter);
        end

    end
end

if mpc.has_h_cnstr

    if h_cnstr.min_limit && h_cnstr.max_limit
        if h_cnstr.use_k0
            row_min = h_cnstr.min_ineqRow_0;
            row_max = h_cnstr.max_ineqRow_0;

            iS_ri_hat_h0 = iS_ri_hat_0(row_max) - iS_ri_hat_0(row_min);
            iS_h0 = iS_0(row_min) + iS_0(row_max);

            ru_0 = ru_0 + Dh_0' * iS_ri_hat_h0;
            R_0 = R_0 + Dh_0' * (iS_h0 .* Dh_0);
        end

        % k = 1,...,N-1
        row_min = h_cnstr.min_ineqRow_k;
        row_max = h_cnstr.max_ineqRow_k;
        for k = 1:N-1

            iS_ri_hat_hk = iS_ri_hat_k(row_max,k) - iS_ri_hat_k(row_min,k);
            iS_hk = iS_k(row_min,k) + iS_k(row_max,k);

            if h_cnstr.use_s && h_cnstr.use_su && h_cnstr.use_u

                rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_hk;
                rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_hk;
                ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_hk;

                iS_Ch_k = iS_hk .* Ch(:,:,k);
                iS_Dsuh_k = iS_hk .* Dsuh(:,:,k);

                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
                Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
                Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_hk .* Dh(:,:,k));
                Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
                Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
            elseif h_cnstr.use_s && h_cnstr.use_su

                rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_hk;
                rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_hk;

                iS_Ch_k = iS_hk .* Ch(:,:,k);
                iS_Dsuh_k = iS_hk .* Dsuh(:,:,k);

                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
                Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
                Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
            elseif h_cnstr.use_s && h_cnstr.use_u

                rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_hk;
                ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_hk;

                iS_Ch_k = iS_hk .* Ch(:,:,k);

                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_hk .* Dh(:,:,k));
                Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
            elseif h_cnstr.use_su && h_cnstr.use_u

                rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_hk;
                ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_hk;

                iS_Dsuh_k = iS_hk .* Dsuh(:,:,k);

                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_hk .* Dh(:,:,k));
                Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
            elseif h_cnstr.use_s

                rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_hk;
                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * (iS_hk .* Ch(:,:,k));
            elseif h_cnstr.use_su

                rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_hk;
                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * (iS_hk .* Dsuh(:,:,k));
            elseif h_cnstr.use_u

                ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_hk;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_hk .* Dh(:,:,k));
            end
        end

        % k = N
        if h_cnstr.use_ter
            row_min = h_cnstr.min_ineqRow_ter;
            row_max = h_cnstr.max_ineqRow_ter;

            iS_ri_hat_hter = iS_ri_hat_ter(row_max) - iS_ri_hat_ter(row_min);
            iS_hter = iS_ter(row_min) + iS_ter(row_max);

            rse_ter(s_col) = rse_ter(s_col) + Ch_ter' * iS_ri_hat_hter;
            Q_ter(s_col,s_col) = Q_ter(s_col,s_col) + Ch_ter' * (iS_hter .* Ch_ter);
        end

    elseif h_cnstr.min_limit

        if h_cnstr.use_k0
            row_min = h_cnstr.min_ineqRow_0;

            ru_0 = ru_0 - Dh_0' * iS_ri_hat_0(row_min);
            R_0 = R_0 + Dh_0' * (iS_0(row_min) .* Dh_0);
        end

        % k = 1,...,N-1
        row_min = h_cnstr.min_ineqRow_k;
        for k = 1:N-1

            if h_cnstr.use_s && h_cnstr.use_su && h_cnstr.use_u

                rse_k(s_col,k) = rse_k(s_col,k) - Ch(:,:,k)' * iS_ri_hat_k(row_min,k);
                rse_k(su_col,k) = rse_k(su_col,k) - Dsuh(:,:,k)' * iS_ri_hat_k(row_min,k);
                ru_k(:,k) = ru_k(:,k) - Dh(:,:,k)' * iS_ri_hat_k(row_min,k);

                iS_Ch_k = iS_k(row_min,k) .* Ch(:,:,k);
                iS_Dsuh_k = iS_k(row_min,k) .* Dsuh(:,:,k);

                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
                Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
                Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(row_min,k) .* Dh(:,:,k));
                Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
                Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
            elseif h_cnstr.use_s && h_cnstr.use_su

                rse_k(s_col,k) = rse_k(s_col,k) - Ch(:,:,k)' * iS_ri_hat_k(row_min,k);
                rse_k(su_col,k) = rse_k(su_col,k) - Dsuh(:,:,k)' * iS_ri_hat_k(row_min,k);

                iS_Ch_k = iS_k(row_min,k) .* Ch(:,:,k);
                iS_Dsuh_k = iS_k(row_min,k) .* Dsuh(:,:,k);

                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
                Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
                Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
            elseif h_cnstr.use_s && h_cnstr.use_u

                rse_k(s_col,k) = rse_k(s_col,k) - Ch(:,:,k)' * iS_ri_hat_k(row_min,k);
                ru_k(:,k) = ru_k(:,k) - Dh(:,:,k)' * iS_ri_hat_k(row_min,k);

                iS_Ch_k = iS_k(row_min,k) .* Ch(:,:,k);

                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(row_min,k) .* Dh(:,:,k));
                Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
            elseif h_cnstr.use_su && h_cnstr.use_u

                rse_k(su_col,k) = rse_k(su_col,k) - Dsuh(:,:,k)' * iS_ri_hat_k(row_min,k);
                ru_k(:,k) = ru_k(:,k) - Dh(:,:,k)' * iS_ri_hat_k(row_min,k);

                iS_Dsuh_k = iS_k(row_min,k) .* Dsuh(:,:,k);

                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(row_min,k) .* Dh(:,:,k));
                Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
            elseif h_cnstr.use_s

                rse_k(s_col,k) = rse_k(s_col,k) - Ch(:,:,k)' * iS_ri_hat_k(row_min,k);
                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * (iS_k(row_min,k) .* Ch(:,:,k));
            elseif h_cnstr.use_su

                rse_k(su_col,k) = rse_k(su_col,k) - Dsuh(:,:,k)' * iS_ri_hat_k(row_min,k);
                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * (iS_k(row_min,k) .* Dsuh(:,:,k));
            elseif h_cnstr.use_u

                ru_k(:,k) = ru_k(:,k) - Dh(:,:,k)' * iS_ri_hat_k(row_min,k);
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(row_min,k) .* Dh(:,:,k));
            end
        end

        % k = N
        if h_cnstr.use_ter
            row_min = h_cnstr.min_ineqRow_ter;

            rse_ter(s_col) = rse_ter(s_col) - Ch_ter' * iS_ri_hat_ter(row_min);
            Q_ter(s_col,s_col) = Q_ter(s_col,s_col) + Ch_ter' * (iS_ter(row_min) .* Ch_ter);
        end

    elseif h_cnstr.max_limit

        if h_cnstr.use_k0
            row_max = h_cnstr.max_ineqRow_0;

            ru_0 = ru_0 + Dh_0' * iS_ri_hat_0(row_max);
            R_0 = R_0 + Dh_0' * (iS_0(row_max) .* Dh_0);
        end

        % k = 1,...,N-1
        row_max = h_cnstr.max_ineqRow_k;
        for k = 1:N-1

            if h_cnstr.use_s && h_cnstr.use_su && h_cnstr.use_u

                rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_k(row_max,k);
                rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_k(row_max,k);
                ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_k(row_max,k);

                iS_Ch_k = iS_k(row_max,k) .* Ch(:,:,k);
                iS_Dsuh_k = iS_k(row_max,k) .* Dsuh(:,:,k);

                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
                Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
                Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(row_max,k) .* Dh(:,:,k));
                Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
                Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
            elseif h_cnstr.use_s && h_cnstr.use_su

                rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_k(row_max,k);
                rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_k(row_max,k);

                iS_Ch_k = iS_k(row_max,k) .* Ch(:,:,k);
                iS_Dsuh_k = iS_k(row_max,k) .* Dsuh(:,:,k);

                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
                Q_k(s_col,su_col,k) = Q_k(s_col,su_col,k) + Ch(:,:,k)' * iS_Dsuh_k;
                Q_k(su_col,s_col,k) = Q_k(su_col,s_col,k) + Dsuh(:,:,k)' * iS_Ch_k;
                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
            elseif h_cnstr.use_s && h_cnstr.use_u

                rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_k(row_max,k);
                ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_k(row_max,k);

                iS_Ch_k = iS_k(row_max,k) .* Ch(:,:,k);

                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * iS_Ch_k;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(row_max,k) .* Dh(:,:,k));
                Y_k(:,s_col,k) = Y_k(:,s_col,k) + Dh(:,:,k)' * iS_Ch_k;
            elseif h_cnstr.use_su && h_cnstr.use_u

                rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_k(row_max,k);
                ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_k(row_max,k);

                iS_Dsuh_k = iS_k(row_max,k) .* Dsuh(:,:,k);

                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * iS_Dsuh_k;
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(row_max,k) .* Dh(:,:,k));
                Y_k(:,su_col,k) = Y_k(:,su_col,k) + Dh(:,:,k)' * iS_Dsuh_k;
            elseif h_cnstr.use_s

                rse_k(s_col,k) = rse_k(s_col,k) + Ch(:,:,k)' * iS_ri_hat_k(row_max,k);
                Q_k(s_col,s_col,k) = Q_k(s_col,s_col,k) + Ch(:,:,k)' * (iS_k(row_max,k) .* Ch(:,:,k));
            elseif h_cnstr.use_su

                rse_k(su_col,k) = rse_k(su_col,k) + Dsuh(:,:,k)' * iS_ri_hat_k(row_max,k);
                Q_k(su_col,su_col,k) = Q_k(su_col,su_col,k) + Dsuh(:,:,k)' * (iS_k(row_max,k) .* Dsuh(:,:,k));
            elseif h_cnstr.use_u

                ru_k(:,k) = ru_k(:,k) + Dh(:,:,k)' * iS_ri_hat_k(row_max,k);
                R_k(:,:,k) = R_k(:,:,k) + Dh(:,:,k)' * (iS_k(row_max,k) .* Dh(:,:,k));
            end
        end

        % k = N
        if h_cnstr.use_ter
            row_max = h_cnstr.max_ineqRow_ter;

            rse_ter(s_col) = rse_ter(s_col) + Ch_ter' * iS_ri_hat_ter(row_max);
            Q_ter(s_col,s_col) = Q_ter(s_col,s_col) + Ch_ter' * (iS_ter(row_max) .* Ch_ter);
        end

    end
end

end
