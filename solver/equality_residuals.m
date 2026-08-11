function mpc = equality_residuals(mpc)

% k = 0
mpc.rp_0(mpc.s_col) = mpc.B(:,:,1)*mpc.u(:,1)-mpc.s(:,1)-mpc.beq_0;

for k = 1:mpc.N-1
    % [A B -I][s u s+]'-beq
    mpc.rp_k(mpc.s_col,k) = mpc.A(:,:,k+1)*mpc.s(:,k)+mpc.B(:,:,k+1)*mpc.u(:,k+1)-mpc.s(:,k+1)...
                            -mpc.beq_k(:,k);
end

if mpc.has_du
    % this is the control input delay equality condition: su+ = u_prev
    % [I -I][u su+]' = 0
    % k = 0
    mpc.rp_0(mpc.su_col) = mpc.u(:,1)-mpc.su(:,1);
    % k = 1...N
    mpc.rp_k(mpc.su_col,:) = mpc.u(:,2:mpc.N)-mpc.su(:,2:mpc.N);
end

%%
% Inequality residuals are assembled directly from materialized signals.
% A lower-bound row is q_min - q + g - E*v and an upper-bound row is
% q - q_max + g - E*v.

% k = 0
mpc.ri_0(:) = mpc.g_0;
if mpc.nv_k(1)
    mpc.ri_0(mpc.v_rows_0) = mpc.ri_0(mpc.v_rows_0) - mpc.v_0;
end

if mpc.has_u_cnstr
    if mpc.u_cnstr.min_limit
        row = mpc.u_cnstr.min_ineqRow_0;
        mpc.ri_0(row) = mpc.ri_0(row) + mpc.u_cnstr.min(:,1) - mpc.u(:,1);
    end
    if mpc.u_cnstr.max_limit
        row = mpc.u_cnstr.max_ineqRow_0;
        mpc.ri_0(row) = mpc.ri_0(row) + mpc.u(:,1) - mpc.u_cnstr.max(:,1);
    end
end

if mpc.has_du_cnstr
    if mpc.du_cnstr.min_limit
        row = mpc.du_cnstr.min_ineqRow_0;
        mpc.ri_0(row) = mpc.ri_0(row) + mpc.du_cnstr.min(:,1) - mpc.du(:,1);
    end
    if mpc.du_cnstr.max_limit
        row = mpc.du_cnstr.max_ineqRow_0;
        mpc.ri_0(row) = mpc.ri_0(row) + mpc.du(:,1) - mpc.du_cnstr.max(:,1);
    end
end

if mpc.has_y_cnstr
    if mpc.y_cnstr.min_limit && mpc.y_use_k0
        row = mpc.y_cnstr.min_ineqRow_0;
        mpc.ri_0(row) = mpc.ri_0(row) + mpc.y_cnstr.min_0 - mpc.y_0;
    end
    if mpc.y_cnstr.max_limit && mpc.y_use_k0
        row = mpc.y_cnstr.max_ineqRow_0;
        mpc.ri_0(row) = mpc.ri_0(row) + mpc.y_0 - mpc.y_cnstr.max_0;
    end
end

if mpc.has_h_cnstr
    if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_k0
        row = mpc.h_cnstr.min_ineqRow_0;
        mpc.ri_0(row) = mpc.ri_0(row) + mpc.h_cnstr.min_0 - mpc.h_0;
    end
    if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_k0
        row = mpc.h_cnstr.max_ineqRow_0;
        mpc.ri_0(row) = mpc.ri_0(row) + mpc.h_0 - mpc.h_cnstr.max_0;
    end
end

if mpc.N > 1
    mpc.ri_k(:,:) = mpc.g_k;
    if mpc.nv_k(2)
        mpc.ri_k(mpc.v_rows_k,:) = mpc.ri_k(mpc.v_rows_k,:) - mpc.v_k;
    end

    if mpc.has_s_cnstr
        if mpc.s_cnstr.min_limit
            row = mpc.s_cnstr.min_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.s_cnstr.min(:,1:mpc.N-1) - mpc.s(:,1:mpc.N-1);
        end
        if mpc.s_cnstr.max_limit
            row = mpc.s_cnstr.max_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.s(:,1:mpc.N-1) - mpc.s_cnstr.max(:,1:mpc.N-1);
        end
    end

    if mpc.has_u_cnstr
        if mpc.u_cnstr.min_limit
            row = mpc.u_cnstr.min_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.u_cnstr.min(:,2:mpc.N) - mpc.u(:,2:mpc.N);
        end
        if mpc.u_cnstr.max_limit
            row = mpc.u_cnstr.max_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.u(:,2:mpc.N) - mpc.u_cnstr.max(:,2:mpc.N);
        end
    end

    if mpc.has_du_cnstr
        if mpc.du_cnstr.min_limit
            row = mpc.du_cnstr.min_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.du_cnstr.min(:,2:mpc.N) - mpc.du(:,2:mpc.N);
        end
        if mpc.du_cnstr.max_limit
            row = mpc.du_cnstr.max_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.du(:,2:mpc.N) - mpc.du_cnstr.max(:,2:mpc.N);
        end
    end

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit
            row = mpc.y_cnstr.min_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.y_cnstr.min(:,1:mpc.N-1) - mpc.y(:,1:mpc.N-1);
        end
        if mpc.y_cnstr.max_limit
            row = mpc.y_cnstr.max_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.y(:,1:mpc.N-1) - mpc.y_cnstr.max(:,1:mpc.N-1);
        end
    end

    if mpc.has_h_cnstr
        if mpc.h_cnstr.min_limit
            row = mpc.h_cnstr.min_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.h_cnstr.min(:,1:mpc.N-1) - mpc.h(:,1:mpc.N-1);
        end
        if mpc.h_cnstr.max_limit
            row = mpc.h_cnstr.max_ineqRow_k;
            mpc.ri_k(row,:) = mpc.ri_k(row,:) + ...
                mpc.h(:,1:mpc.N-1) - mpc.h_cnstr.max(:,1:mpc.N-1);
        end
    end
end

if mpc.ng_k(3)
    mpc.ri_ter(:) = mpc.g_ter - mpc.v_ter;

    if mpc.has_s_cnstr
        if mpc.s_cnstr.min_limit
            row = mpc.s_cnstr.min_ineqRow_ter;
            mpc.ri_ter(row) = mpc.ri_ter(row) + mpc.s_cnstr.min(:,mpc.N) - mpc.s(:,mpc.N);
        end
        if mpc.s_cnstr.max_limit
            row = mpc.s_cnstr.max_ineqRow_ter;
            mpc.ri_ter(row) = mpc.ri_ter(row) + mpc.s(:,mpc.N) - mpc.s_cnstr.max(:,mpc.N);
        end
    end

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit && mpc.y_use_ter
            row = mpc.y_cnstr.min_ineqRow_ter;
            mpc.ri_ter(row) = mpc.ri_ter(row) + mpc.y_cnstr.min_ter - mpc.y_ter;
        end
        if mpc.y_cnstr.max_limit && mpc.y_use_ter
            row = mpc.y_cnstr.max_ineqRow_ter;
            mpc.ri_ter(row) = mpc.ri_ter(row) + mpc.y_ter - mpc.y_cnstr.max_ter;
        end
    end

    if mpc.has_h_cnstr
        if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_ter
            row = mpc.h_cnstr.min_ineqRow_ter;
            mpc.ri_ter(row) = mpc.ri_ter(row) + mpc.h_cnstr.min_ter - mpc.h_ter;
        end
        if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_ter
            row = mpc.h_cnstr.max_ineqRow_ter;
            mpc.ri_ter(row) = mpc.ri_ter(row) + mpc.h_ter - mpc.h_cnstr.max_ter;
        end
    end
end

end
