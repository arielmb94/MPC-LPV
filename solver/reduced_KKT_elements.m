function mpc = reduced_KKT_elements(mpc)
%% ri_hat = ri - g^2*(-1/g) + v^2*(t*qv-1/v)
%  ri_hat = ri + g + t*qv*v^2 - v


if mpc.ng_k(1)
    mpc.ri_hat_0 = mpc.ri_0+mpc.g_0;
    mpc.g2_0 = mpc.g_0.^2;
    if mpc.nv_k(1)
        mpc.v2_0 = mpc.v_0.^2;
        % v^2*rv = t*qv*v^2 - v
        mpc.rv_v2_0 = mpc.t*mpc.grad_qv_0.*mpc.v2_0-mpc.v_0;
        mpc.ri_hat_0 = mpc.ri_hat_0(mpc.v_rows_0) + mpc.rv_v2_0;
    end
end

if mpc.ng_k(2)
    mpc.ri_hat_k = mpc.ri_k+mpc.g_k;
    mpc.g2_k = mpc.g_k.^2;
    if mpc.nv_k(2)
        mpc.v2_k = mpc.v_k.^2;
        % v^2*rv = t*qv*v^2 - v
        mpc.rv_v2_k = mpc.t*mpc.grad_qv_k.*mpc.v2_k-mpc.v_k;
        mpc.ri_hat_k(mpc.v_rows_k,:) = mpc.ri_hat_k(mpc.v_rows_k,:) + mpc.rv_v2_k;
    end
end

if mpc.ng_k(3)
    mpc.g2_ter = mpc.g_ter.^2;
    mpc.v2_ter = mpc.v_ter.^2;
    % v^2*rv = t*qv*v^2 - v
    mpc.rv_v2_ter = mpc.t*mpc.grad_qv_ter.*mpc.v2_ter-mpc.v_ter;
    mpc.ri_hat_ter = mpc.ri_ter+mpc.g_ter + mpc.rv_v2_ter;
end

%% S = g^2 + v^2

if mpc.ng_k(1)
    mpc.iS_0 = mpc.g2_0;
    if mpc.nv_k(1), mpc.iS_0(mpc.v_rows_0) = mpc.iS_0(mpc.v_rows_0) + mpc.v2_0; end
    mpc.iS_0  = 1./mpc.iS_0;
    mpc.iS_ri_hat_0 = mpc.iS_0.*mpc.ri_hat_0;
end

if mpc.ng_k(2)
    mpc.iS_k = mpc.g2_k;
    if mpc.nv_k(2), mpc.iS_k(mpc.v_rows_k,:) = mpc.iS_k(mpc.v_rows_k,:) + mpc.v2_k; end
    mpc.iS_k  = 1./mpc.iS_k;
    mpc.iS_ri_hat_k = mpc.iS_k.*mpc.ri_hat_k;
end

if mpc.ng_k(3)
    mpc.iS_ter = mpc.g2_ter;
    if mpc.nv_k(3), mpc.iS_ter = mpc.iS_ter + mpc.v2_ter; end
    mpc.iS_ter  = 1./mpc.iS_ter;
    mpc.iS_ri_hat_ter = mpc.iS_ter.*mpc.ri_hat_ter;
end

%% rx_hat

% Cache the arrays used by the constraint-major updates.
if mpc.ng_k(1)
    iS_0 = mpc.iS_0;
    iS_ri_hat_0 = mpc.iS_ri_hat_0;
end
if mpc.ng_k(2)
    iS_k = mpc.iS_k;
    iS_ri_hat_k = mpc.iS_ri_hat_k;
end
if mpc.ng_k(3)
    iS_ter = mpc.iS_ter;
    iS_ri_hat_ter = mpc.iS_ri_hat_ter;
end
s_col = mpc.s_col;
su_col = mpc.su_col;

% k = 0
%mpc.iS_ri_hat_0 = iS_0.*ri_hat_0;

ru_hat_0 = mpc.ru_0;

R_0 = mpc.t*mpc.H_f0_0;

% k = 1,...,N-1
%mpc.iS_ri_hat_k = iS_k.*ri_hat_k;

rse_hat_k = mpc.rse_k;
ru_hat_k = mpc.ru_k;

R_k = mpc.t*mpc.H_f0_k(mpc.u_col,mpc.u_col,:);
Q_k = mpc.t*mpc.H_f0_k(mpc.se_col,mpc.se_col,:);
Y_k = mpc.t*mpc.H_f0_k(mpc.u_col,mpc.se_col,:);

% k = N
%mpc.iS_ri_hat_ter = iS_ter.*ri_hat_ter;

rse_hat_ter = mpc.rse_ter;
Q_ter = mpc.t*mpc.H_f0_ter;

if mpc.has_u_cnstr
    if mpc.u_cnstr.min_limit && mpc.u_cnstr.max_limit

        row_min = mpc.u_cnstr.min_ineqRow_0;
        row_max = mpc.u_cnstr.max_ineqRow_0;

        % k = 0
        for i = 1:mpc.nu
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            iS_min_i = iS_0(row_min_i);
            iS_max_i = iS_0(row_max_i);
            iS_ri_hat_min_i = iS_ri_hat_0(row_min_i);
            iS_ri_hat_max_i = iS_ri_hat_0(row_max_i);

            ru_hat_0(i) = ru_hat_0(i) + ...
                iS_ri_hat_max_i - iS_ri_hat_min_i;
            R_0(i,i) = R_0(i,i) + iS_min_i + iS_max_i;

        end

        % k = 1,...,N-1
        row_min = mpc.u_cnstr.min_ineqRow_k;
        row_max = mpc.u_cnstr.max_ineqRow_k;

        for k = 1:mpc.N-1
            for i = 1:mpc.nu
                row_min_i = row_min(i);
                row_max_i = row_max(i);

                iS_min_i = iS_k(row_min_i,k);
                iS_max_i = iS_k(row_max_i,k);
                iS_ri_hat_min_i = iS_ri_hat_k(row_min_i,k);
                iS_ri_hat_max_i = iS_ri_hat_k(row_max_i,k);

                ru_hat_k(i,k) = ru_hat_k(i,k) + ...
                    iS_ri_hat_max_i - iS_ri_hat_min_i;
                R_k(i,i,k) = R_k(i,i,k) + iS_min_i + iS_max_i;

            end
        end
    elseif mpc.u_cnstr.min_limit
        
        % k = 0
        row = mpc.u_cnstr.min_ineqRow_0;

        for i = 1:mpc.nu
            row_i = row(i);

            iS_i = iS_0(row_i);
            iS_ri_hat_i = iS_ri_hat_0(row_i);

            % A_i = -I <-- u
            ru_hat_0(i) = ru_hat_0(i) - iS_ri_hat_i;
            R_0(i,i) =  R_0(i,i) + iS_i;

        end

        % k = 1,...,N-1
        row = mpc.u_cnstr.min_ineqRow_k;

        for k = 1:mpc.N-1
            for i = 1:mpc.nu

                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = -I <-- u
                ru_hat_k(i,k) = ru_hat_k(i,k) - iS_ri_hat_i;

                R_k(i,i,k) =  R_k(i,i,k) + iS_i;

            end
        end
    elseif mpc.u_cnstr.max_limit
        
        % k = 0
        row = mpc.u_cnstr.max_ineqRow_0;

        for i = 1:mpc.nu
            row_i = row(i);

            iS_i = iS_0(row_i);
            iS_ri_hat_i = iS_ri_hat_0(row_i);

            % A_i = I <-- u
            ru_hat_0(i) = ru_hat_0(i) + iS_ri_hat_i;
            R_0(i,i) =  R_0(i,i) + iS_i;

        end

        % k = 1,...,N-1
        row = mpc.u_cnstr.max_ineqRow_k;

        for k = 1:mpc.N-1
            for i = 1:mpc.nu

                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = I <-- u
                ru_hat_k(i,k) = ru_hat_k(i,k) + iS_ri_hat_i;

                R_k(i,i,k) =  R_k(i,i,k) + iS_i;

            end
        end
    end
end

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit && mpc.du_cnstr.max_limit

        row_min = mpc.du_cnstr.min_ineqRow_0;
        row_max = mpc.du_cnstr.max_ineqRow_0;

        % k = 0
        for i = 1:mpc.nu
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            iS_min_i = iS_0(row_min_i);
            iS_max_i = iS_0(row_max_i);
            iS_ri_hat_min_i = iS_ri_hat_0(row_min_i);
            iS_ri_hat_max_i = iS_ri_hat_0(row_max_i);

            ru_hat_0(i) = ru_hat_0(i) + ...
                iS_ri_hat_max_i - iS_ri_hat_min_i;
            R_0(i,i) = R_0(i,i) + iS_min_i + iS_max_i;

        end

        % k = 1,...,N-1
        row_min = mpc.du_cnstr.min_ineqRow_k;
        row_max = mpc.du_cnstr.max_ineqRow_k;

        for k = 1:mpc.N-1
            for i = 1:mpc.nu
                su_i = su_col(i);
                row_min_i = row_min(i);
                row_max_i = row_max(i);

                iS_min_i = iS_k(row_min_i,k);
                iS_max_i = iS_k(row_max_i,k);
                iS_ri_hat_min_i = iS_ri_hat_k(row_min_i,k);
                iS_ri_hat_max_i = iS_ri_hat_k(row_max_i,k);

                iS_ri_hat_delta_i = iS_ri_hat_max_i - iS_ri_hat_min_i;
                iS_sum_i = iS_min_i + iS_max_i;

                rse_hat_k(su_i,k) = rse_hat_k(su_i,k) - ...
                    iS_ri_hat_delta_i;
                ru_hat_k(i,k) = ru_hat_k(i,k) + iS_ri_hat_delta_i;

                Q_k(su_i,su_i,k) = Q_k(su_i,su_i,k) + iS_sum_i;
                R_k(i,i,k) = R_k(i,i,k) + iS_sum_i;
                Y_k(i,su_i,k) = Y_k(i,su_i,k) - iS_sum_i;

            end
        end
    elseif mpc.du_cnstr.min_limit
        
        % k = 0
        row = mpc.du_cnstr.min_ineqRow_0;

        for i = 1:mpc.nu
            row_i = row(i);

            iS_i = iS_0(row_i);
            iS_ri_hat_i = iS_ri_hat_0(row_i);

            % A_i = -I <-- u
            ru_hat_0(i) = ru_hat_0(i) - iS_ri_hat_i;
            R_0(i,i) =  R_0(i,i) + iS_i;

        end

        % k = 1,...,N-1
        row = mpc.du_cnstr.min_ineqRow_k;

        for k = 1:mpc.N-1
            for i = 1:mpc.nu

                su_i = su_col(i);
                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = [I -I] <-- [su u]^T
                rse_hat_k(su_i,k) = rse_hat_k(su_i,k) + iS_ri_hat_i;
                ru_hat_k(i,k) = ru_hat_k(i,k) - iS_ri_hat_i;

                Q_k(su_i,su_i,k) = Q_k(su_i,su_i,k) + iS_i;
                R_k(i,i,k) =  R_k(i,i,k) + iS_i;
                Y_k(i,su_i,k) = Y_k(i,su_i,k) - iS_i;

            end
        end
    elseif mpc.du_cnstr.max_limit
        
        % k = 0
        row = mpc.du_cnstr.max_ineqRow_0;

        for i = 1:mpc.nu
            row_i = row(i);

            iS_i = iS_0(row_i);
            iS_ri_hat_i = iS_ri_hat_0(row_i);

            % A_i = I <-- u
            ru_hat_0(i) = ru_hat_0(i) + iS_ri_hat_i;
            R_0(i,i) =  R_0(i,i) + iS_i;

        end

        % k = 1,...,N-1
        row = mpc.du_cnstr.max_ineqRow_k;

        for k = 1:mpc.N-1
            for i = 1:mpc.nu

                su_i = su_col(i);
                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = [-I I] <-- [su u]^T
                rse_hat_k(su_i,k) = rse_hat_k(su_i,k) - iS_ri_hat_i;
                ru_hat_k(i,k) = ru_hat_k(i,k) + iS_ri_hat_i;

                Q_k(su_i,su_i,k) = Q_k(su_i,su_i,k) + iS_i;
                R_k(i,i,k) =  R_k(i,i,k) + iS_i;
                Y_k(i,su_i,k) = Y_k(i,su_i,k) - iS_i;

            end
        end
    end
end

if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit && mpc.s_cnstr.max_limit

        row_min = mpc.s_cnstr.min_ineqRow_k;
        row_max = mpc.s_cnstr.max_ineqRow_k;

        % k = 1,...,N-1
        for k = 1:mpc.N-1
            for i = 1:mpc.nx
                s_i = s_col(i);
                row_min_i = row_min(i);
                row_max_i = row_max(i);

                iS_min_i = iS_k(row_min_i,k);
                iS_max_i = iS_k(row_max_i,k);
                iS_ri_hat_min_i = iS_ri_hat_k(row_min_i,k);
                iS_ri_hat_max_i = iS_ri_hat_k(row_max_i,k);

                rse_hat_k(s_i,k) = rse_hat_k(s_i,k) + ...
                    iS_ri_hat_max_i - iS_ri_hat_min_i;
                Q_k(s_i,s_i,k) = Q_k(s_i,s_i,k) + ...
                    iS_min_i + iS_max_i;

            end
        end

        % k = N
        row_min = mpc.s_cnstr.min_ineqRow_ter;
        row_max = mpc.s_cnstr.max_ineqRow_ter;

        for i = 1:mpc.nx
            row_min_i = row_min(i);
            row_max_i = row_max(i);

            iS_min_i = iS_ter(row_min_i);
            iS_max_i = iS_ter(row_max_i);
            iS_ri_hat_min_i = iS_ri_hat_ter(row_min_i);
            iS_ri_hat_max_i = iS_ri_hat_ter(row_max_i);

            rse_hat_ter(i) = rse_hat_ter(i) + ...
                iS_ri_hat_max_i - iS_ri_hat_min_i;
            Q_ter(i,i) = Q_ter(i,i) + iS_min_i + iS_max_i;

        end
    elseif mpc.s_cnstr.min_limit

        % k = 1,...,N-1
        row = mpc.s_cnstr.min_ineqRow_k;

        for k = 1:mpc.N-1
            for i = 1:mpc.nx

                s_i = s_col(i);
                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = -I <-- s
                rse_hat_k(s_i,k) = rse_hat_k(s_i,k) - iS_ri_hat_i;

                Q_k(s_i,s_i,k) = Q_k(s_i,s_i,k) + iS_i;

            end
        end

        % k = N
        row = mpc.s_cnstr.min_ineqRow_ter;

        for i = 1:mpc.nx
            row_i = row(i);

            iS_i = iS_ter(row_i);
            iS_ri_hat_i = iS_ri_hat_ter(row_i);

            % A_i = -I <-- s
            rse_hat_ter(i) = rse_hat_ter(i) - iS_ri_hat_i;
            Q_ter(i,i) =  Q_ter(i,i) + iS_i;

        end
    elseif mpc.s_cnstr.max_limit

        % k = 1,...,N-1
        row = mpc.s_cnstr.max_ineqRow_k;

        for k = 1:mpc.N-1
            for i = 1:mpc.nx

                s_i = s_col(i);
                row_i = row(i);

                iS_i = iS_k(row_i,k);
                iS_ri_hat_i = iS_ri_hat_k(row_i,k);

                % A_i = I <-- s
                rse_hat_k(s_i,k) = rse_hat_k(s_i,k) + iS_ri_hat_i;

                Q_k(s_i,s_i,k) = Q_k(s_i,s_i,k) + iS_i;

            end
        end

        % k = N
        row = mpc.s_cnstr.max_ineqRow_ter;

        for i = 1:mpc.nx
            row_i = row(i);

            iS_i = iS_ter(row_i);
            iS_ri_hat_i = iS_ri_hat_ter(row_i);

            % A_i = I <-- s
            rse_hat_ter(i) = rse_hat_ter(i) + iS_ri_hat_i;
            Q_ter(i,i) =  Q_ter(i,i) + iS_i;

        end
    end
end


mpc.ru_hat_0 = ru_hat_0;
mpc.R_0 = R_0;
mpc.rse_hat_k = rse_hat_k;
mpc.ru_hat_k = ru_hat_k;
mpc.R_k = R_k;
mpc.Q_k = Q_k;
mpc.Y_k = Y_k;
mpc.rse_hat_ter = rse_hat_ter;
mpc.Q_ter = Q_ter;
% if mpc.has_y_cnstr && mpc.y_use_k0
%     if mpc.y_cnstr.min_limit
%         row = mpc.y_cnstr.min_ineqRow_0;
%         mpc.ru_hat_0 = mpc.ru_hat_0 - mpc.D_0'*mpc.iS_ri_hat_0(row);
%     end
%     if mpc.y_cnstr.max_limit
%         row = mpc.y_cnstr.max_ineqRow_0;
%         mpc.ru_hat_0 = mpc.ru_hat_0 + mpc.D_0'*mpc.iS_ri_hat_0(row);
%     end
% end
% 
% if mpc.has_h_cnstr 
%     if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_k0
%         row = mpc.h_cnstr.min_ineqRow_0;
%         mpc.ru_hat_0 = mpc.ru_hat_0 - mpc.Dh_0'*mpc.iS_ri_hat_0(row);
%     end
%     if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_k0
%         row = mpc.h_cnstr.max_ineqRow_0;
%         mpc.ru_hat_0 = mpc.ru_hat_0 + mpc.Dh_0'*mpc.iS_ri_hat_0(row);
%     end
% end



% 
% for k = 1:mpc.N-1
% 
%     mpc.rse_hat_k(:,k) = mpc.rse_k(:,k) + mpc.Ai_k(:,mpc.se_col,k)'*mpc.iS_ri_hat_k(:,k);
%     mpc.ru_hat_k(:,k) = mpc.ru_k(:,k) + mpc.Ai_k(:,mpc.u_col,k)'*mpc.iS_ri_hat_k(:,k);
% end

%mpc.rse_hat_ter = mpc.rse_ter + mpc.Ai_ter'*mpc.iS_ri_hat_ter;

%% H = Hess(f0) + Ai'*(S^-1)*Ai

% mpc.iS_Ai_0 = mpc.iS_0.*mpc.Ai_0;
% %mpc.R_0(:,:) = mpc.t*mpc.H_f0_0 + mpc.Ai_0'*mpc.iS_Ai_0;
% 
% for k = 1:mpc.N-1
% 
%     mpc.iS_Ai_k(:,:,k) = mpc.iS_k(:,k).*mpc.Ai_k(:,:,k);
% 
% %     mpc.H_k(:,:,k) = mpc.t*mpc.H_f0_k(:,:,k) + ...
% %                         mpc.Ai_k(:,:,k)'*mpc.iS_Ai_k(:,:,k);
% % 
% %     mpc.Q_k(:,:,k) = mpc.H_k(mpc.se_col,mpc.se_col,k);
% %     mpc.R_k(:,:,k) = mpc.H_k(mpc.u_col,mpc.u_col,k);
% %     mpc.Y_k(:,:,k) = mpc.H_k(mpc.u_col,mpc.se_col,k);
% end
% 
% %mpc.Q_ter(:,:) = mpc.t*mpc.H_f0_ter;
% if mpc.ng_k(3)
%     mpc.iS_Ai_ter = mpc.iS_ter.*mpc.Ai_ter; 
%     %mpc.Q_ter(:,:) = mpc.Q_ter(:,:) + mpc.Ai_ter'*mpc.iS_Ai_ter;
% end

end
