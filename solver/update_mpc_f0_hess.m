function mpc = update_mpc_f0_hess(mpc)

mpc.recompute_cost_hess = 0;

% Reset Hessian of Cost Function to 0 
mpc.H_f0_0(:,:) = 0;
mpc.H_f0_k(:,:,:) = 0;
mpc.H_f0_ter(:,:) = 0;

% k = 0

if mpc.tracking_cost && mpc.y_use_k0
    mpc.H_f0_0 = mpc.H_f0_0 + mpc.H_ErrCost_0;
end
if mpc.quad_control_cost
    mpc.H_f0_0 = mpc.H_f0_0 + mpc.Ru;
end
if mpc.controlrate_cost
    mpc.H_f0_0 = mpc.H_f0_0 + mpc.Rdu;
end
if mpc.quad_custom_cost && mpc.z_use_k0
    mpc.H_f0_0 = mpc.H_f0_0 + mpc.H_CustomCost_0;
end

% k = 1 to N-1

if mpc.tracking_cost
    index = mpc.tracking_cost_index_k;
    mpc.H_f0_k(index, index, :) = mpc.H_f0_k(index, index, :) + mpc.H_ErrCost_k;
end

if mpc.quad_control_cost
    mpc.H_f0_k(mpc.u_col, mpc.u_col, :) = mpc.H_f0_k(mpc.u_col, mpc.u_col, :) + mpc.Ru;
end

if mpc.controlrate_cost
    index = mpc.du_index_k;
    mpc.H_f0_k(index, index, :) = mpc.H_f0_k(index, index, :) + mpc.H_RateCtrl_k;
end

if mpc.quad_custom_cost
    index = mpc.custom_cost_index_k;
    mpc.H_f0_k(index, index, :) = mpc.H_f0_k(index, index, :) + mpc.H_CustomCost_k;
end

% k = N

if mpc.ter_ingredients
    mpc.H_f0_ter(mpc.s_col, mpc.s_col) = mpc.H_f0_ter(mpc.s_col, mpc.s_col) + mpc.P2;
end

if mpc.tracking_cost && mpc.y_use_ter
    mpc.H_f0_ter(mpc.s_col, mpc.s_col) = mpc.H_f0_ter(mpc.s_col, mpc.s_col) + mpc.H_ErrCost_ter;
end

if mpc.quad_custom_cost && mpc.z_use_ter
    mpc.H_f0_ter(mpc.s_col, mpc.s_col) = mpc.H_f0_ter(mpc.s_col, mpc.s_col) + mpc.H_CustomCost_ter;
end

end