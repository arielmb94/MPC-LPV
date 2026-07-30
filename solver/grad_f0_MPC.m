% gradients exclusively on the cost function on s, su and u (no slacks)
function mpc = grad_f0_MPC(mpc)

mpc.grad_f0_0(:) = 0;
mpc.grad_f0_k(:,:) = 0;
mpc.grad_f0_ter(:) = 0;

% k=0
if mpc.quad_control_cost
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.Ru*mpc.u(:,1);
end
if mpc.lin_control_cost
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.ru;
end

if mpc.controlrate_cost
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.Rdu*mpc.du(:,1);
end

if mpc.tracking_cost && mpc.y_use_k0
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.gradErr_Qe_0*mpc.err_0;
end

if mpc.quad_custom_cost && mpc.z_use_k0
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.gradz_Qz_0*mpc.z_0;
end
if mpc.lin_custom_cost && mpc.z_use_k0
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.gradz_qz_0;
end

for k = 1:mpc.N-1

    if mpc.tracking_cost
        index = mpc.tracking_cost_index_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradErr_Qe_k(:,:,k)*mpc.err(:,k);
    end

    if mpc.quad_control_cost
        index = mpc.u_col;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.Ru*mpc.u(:,k+1);
    end
    if mpc.lin_control_cost
        index = mpc.u_col;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.ru;
    end

    if mpc.controlrate_cost
        index = mpc.du_index_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradRateCtrl_Rdu_k(:,:,k)*mpc.du(:,k+1);
    end

    if mpc.quad_custom_cost || mpc.lin_custom_cost
        index = mpc.custom_cost_index_k;

        if mpc.quad_custom_cost
            mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradz_Qz_k(:,:,k)*mpc.z(:,k);
        end
        if mpc.lin_custom_cost
            mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradz_qz_k(:,k);
        end
    end
end

if mpc.ter_ingredients
    mpc.grad_f0_ter(mpc.s_col) = mpc.grad_f0_ter(mpc.s_col) - mpc.P2*(mpc.xN_ref-mpc.s_ter);
end

if mpc.tracking_cost && mpc.y_use_ter
    mpc.grad_f0_ter(mpc.s_col) = mpc.grad_f0_ter(mpc.s_col) + mpc.gradErr_Qe_ter*mpc.err_ter;
end

if mpc.quad_custom_cost && mpc.z_use_ter
    mpc.grad_f0_ter(mpc.s_col) = mpc.grad_f0_ter(mpc.s_col) + mpc.gradz_Qz_ter*mpc.z_ter;
end
if mpc.lin_custom_cost && mpc.z_use_ter
    mpc.grad_f0_ter(mpc.s_col) = mpc.grad_f0_ter(mpc.s_col) + mpc.gradz_qz_ter;
end


mpc.ru_0(:) = mpc.t*mpc.grad_f0_0;

mpc.rse_k(:,:) = mpc.t*mpc.grad_f0_k(mpc.se_col,:);
mpc.ru_k(:,:) = mpc.t*mpc.grad_f0_k(mpc.u_col,:);

mpc.rse_ter(:) = mpc.t*mpc.grad_f0_ter;

end