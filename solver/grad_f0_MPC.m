% gradients exclusively on the cost function on s, su and u (no slacks)
function mpc = grad_f0_MPC(mpc)

mpc.grad_f0_0(:) = 0;
mpc.grad_f0_k(:,:) = 0;
mpc.grad_f0_ter(:) = 0;

% k=0
if mpc.quad_control_cost
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.gradCtlrRu_0*mpc.u(:,1);
end
if mpc.lin_control_cost
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.gradCtlrru_0;
end

if mpc.diffcontrol_cost
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.gradDiffCtlrR_0*mpc.du(:,1);
end

if mpc.tracking_cost && mpc.y_use_k0
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.gradErrQe_0*mpc.err_0;
end

if mpc.quad_custom_cost && mpc.z_use_k0
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.gradPerfQz_0*mpc.z_0;
end
if mpc.lin_custom_cost && mpc.z_use_k0
    mpc.grad_f0_0(:) = mpc.grad_f0_0 + mpc.gradPerfqz_0;
end

for k = 1:mpc.N-1

    if mpc.tracking_cost
        index = mpc.tracking_cost_index_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradErrQe_k(:,:,k)*mpc.err(:,k);
    end

    if mpc.quad_control_cost
        index = mpc.nse+1:mpc.nvar_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradCtlrRu_k(:,:,k)*mpc.u(:,k+1);
    end
    if mpc.lin_control_cost
        index = mpc.nse+1:mpc.nvar_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradCtlrru_k(:,k);
    end

    if mpc.diffcontrol_cost
        index = mpc.du_index_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradDiffCtlrR_k(:,:,k)*mpc.du(:,k+1);
    end

    if mpc.quad_custom_cost || mpc.lin_custom_cost
        index = mpc.custom_cost_index_k;

    if mpc.quad_custom_cost
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradPerfQz_k(:,:,k)*mpc.z(:,k);
    end
    if mpc.lin_custom_cost
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradPerfqz_k(:,k);
    end
    end
end

if mpc.ter_ingredients
    mpc.grad_f0_ter(mpc.s_col) = mpc.grad_f0_ter(mpc.s_col) - mpc.P2*(mpc.xN_ref-mpc.s_ter);
end

if mpc.tracking_cost && mpc.y_use_ter
    mpc.grad_f0_ter(mpc.s_col) = mpc.grad_f0_ter(mpc.s_col) + mpc.gradErrQe_ter*mpc.err_ter;
end

if mpc.quad_custom_cost && mpc.z_use_ter
    mpc.grad_f0_ter(mpc.s_col) = mpc.grad_f0_ter(mpc.s_col) + mpc.gradPerfQz_ter*mpc.z_ter;
end
if mpc.lin_custom_cost && mpc.z_use_ter
    mpc.grad_f0_ter(mpc.s_col) = mpc.grad_f0_ter(mpc.s_col) + mpc.gradPerfqz_ter;
end


mpc.ru_0(:) = mpc.t*mpc.grad_f0_0;

mpc.rse_k(:,:) = mpc.t*mpc.grad_f0_k(mpc.se_col,:);
mpc.ru_k(:,:) = mpc.t*mpc.grad_f0_k(mpc.u_col,:);

mpc.rse_ter(:) = mpc.t*mpc.grad_f0_ter;

end