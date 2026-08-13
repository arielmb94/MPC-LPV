% UPDATE_MPC_OUTPUT_VECTOR Update the tracking-output model.
%
%   mpc = UPDATE_MPC_OUTPUT_VECTOR(mpc, C, D, Dd) updates selected
%   coefficients in
%
%       y_k = C_k*s_k + D_k*u_k + Dd_k*d_k.
%
%   Use [] to leave a coefficient unchanged. The output must first be
%   defined with INIT_MPC_DYNAMICS or INIT_MPC_OUTPUT. This updater cannot
%   change the number of outputs or add a coefficient term that was not
%   enabled during initialization.
%
%   C, D, and Dd may be constant matrices or contain L horizon stages in
%   their third dimension. If L < N, the last supplied stage is reused for
%   the remaining stages. The input d_k is the same fixed known input used
%   by the dynamics.
% 
% 1Existing tracking costs and output constraints use the
%   updated output model.
%
%   Inputs:
%     mpc     - Built CHRONOS MPC structure.
%     C       - Optional state coefficient, size ny-by-nx or ny-by-nx-by-L.
%     D       - Optional control coefficient, size ny-by-nu or
%               ny-by-nu-by-L.
%     Dd      - Optional fixed-known-input coefficient, size ny-by-nd or
%               ny-by-nd-by-L.
%
%   Output:
%     mpc     - Updated CHRONOS MPC structure.
%
%   Example - update only the state and control coefficients:
%
%       mpc = update_mpc_output_vector(mpc, C, D, []);
function mpc = update_mpc_output_vector(mpc,C,D,Dd)

mpc.update_tracking = mpc.tracking_cost;
update_grad = 0;

if ~isempty(C)
    update_grad = 1;

    len_C = size(C,3);
    if len_C < mpc.N
        mpc.C(:,:,:) = fill_mat(mpc.C, C, 1);
        if mpc.y_use_ter, mpc.C_ter(:,:) = mpc.C(mpc.y_rows_ter,:,mpc.N-1); end
    else
        mpc.C(:,:,:) = C(:,:,1:mpc.N-1);
        if mpc.y_use_ter, mpc.C_ter(:,:) = C(mpc.y_rows_ter,:,mpc.N); end
    end
    if mpc.y_use_k0, mpc.C_0(:,:) = mpc.C(mpc.y_rows_k0,:,1); end

    if mpc.tracking_cost && mpc.y_use_ter
        mpc.grad_err_ter(:,:) = -mpc.C_ter';
    end

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit
            for k = 1:mpc.N-1
                mpc.Ai_k(mpc.y_cnstr.min_ineqRow_k,mpc.s_col,k) = -mpc.C(:,:,k);
            end
            if mpc.y_cnstr.use_ter
                mpc.Ai_ter(mpc.y_cnstr.min_ineqRow_ter,mpc.s_col) = -mpc.C_ter;
            end
        end
        if mpc.y_cnstr.max_limit
            for k = 1:mpc.N-1
                mpc.Ai_k(mpc.y_cnstr.max_ineqRow_k,mpc.s_col,k) = mpc.C(:,:,k);
            end
            if mpc.y_cnstr.use_ter
                mpc.Ai_ter(mpc.y_cnstr.max_ineqRow_ter,mpc.s_col) = mpc.C_ter;
            end
        end
    end

end

if ~isempty(D)
    update_grad = 1;

    len_D = size(D,3);
    if len_D < mpc.N-1
        mpc.D(:,:,:) = fill_mat(mpc.D, D, 1);
    else
        mpc.D(:,:,:) = D(:,:,1:mpc.N-1);
    end
    if mpc.y_use_k0
        mpc.D_0(:,:) = mpc.D(mpc.y_rows_k0,:,1);
        if mpc.tracking_cost
            mpc.grad_err_0(:,:) = -mpc.D_0';
        end
    end

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit
            if mpc.y_cnstr.use_k0
                mpc.Ai_0(mpc.y_cnstr.min_ineqRow_0,:) = -mpc.D_0;
            end
            for k = 1:mpc.N-1
                mpc.Ai_k(mpc.y_cnstr.min_ineqRow_k,mpc.u_col,k) = -mpc.D(:,:,k);
            end
        end
        if mpc.y_cnstr.max_limit
            if mpc.y_cnstr.use_k0
                mpc.Ai_0(mpc.y_cnstr.max_ineqRow_0,:) = mpc.D_0;
            end
            for k = 1:mpc.N-1
                mpc.Ai_k(mpc.y_cnstr.max_ineqRow_k,mpc.u_col,k) = mpc.D(:,:,k);
            end
        end
    end
end

if ~isempty(Dd)   

    len_Dd = size(Dd,3);
    if len_Dd < mpc.N-1
        mpc.Dd(:,:,:) = fill_mat(mpc.Dd, Dd, 1);
    else
        mpc.Dd(:,:,:) = Dd(:,:,1:mpc.N-1);
    end

    if mpc.y_use_k0, mpc.Dd_0(:,:) = mpc.Dd(mpc.y_rows_k0,:,1); end
end

if update_grad && mpc.tracking_cost
    for k = 1:mpc.N-1
        if mpc.y_use_s && mpc.y_use_u
            mpc.grad_err(:,:,k) = [-mpc.C(:,:,k)'; -mpc.D(:,:,k)'];
        elseif mpc.y_use_s
            mpc.grad_err(:,:,k) = -mpc.C(:,:,k)';
        elseif mpc.y_use_u
            mpc.grad_err(:,:,k) = -mpc.D(:,:,k)';
        end
    end
end

end
