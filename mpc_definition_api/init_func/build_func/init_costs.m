function mpc = init_costs(mpc)

mpc.H_f0_0 = zeros(mpc.nu);
mpc.H_f0_k = zeros(mpc.nvar_k,mpc.nvar_k,mpc.N-1);
mpc.H_f0_ter = zeros(mpc.nse);

mpc.H_k = zeros(mpc.nvar_k,mpc.nvar_k,mpc.N-1);
mpc.R_0 = zeros(mpc.nu);
mpc.Q_k = zeros(mpc.nse,mpc.nse,mpc.N-1);
mpc.R_k = zeros(mpc.nu,mpc.nu,mpc.N-1);
mpc.Y_k = zeros(mpc.nu,mpc.nse,mpc.N-1);
mpc.Q_ter = zeros(mpc.nse);

mpc.rx_ineq_0 = zeros(mpc.nu,1);
mpc.rx_ineq_k = zeros(mpc.nse+mpc.nu,mpc.N-1);
mpc.rx_ineq_ter = zeros(mpc.nse,1);

mpc.grad_f0_0 = zeros(mpc.nu,1);
mpc.grad_f0_k = zeros(mpc.nvar_k,mpc.N-1);
mpc.grad_f0_ter = zeros(mpc.nse,1);

if mpc.tracking_cost
    mpc = genTrackingCost(mpc);   
end

if mpc.quad_control_cost || mpc.lin_control_cost
    mpc = genControlCost(mpc);
end

if mpc.controlrate_cost
    mpc = genDiffControlCost(mpc);
end

if mpc.quad_custom_cost || mpc.lin_custom_cost
    mpc = genCustomCost(mpc);
end

if any(mpc.nv_k)
   mpc = genSoftSlacksCost(mpc);  
end

if mpc.ter_ingredients
    mpc.P2 = 2*mpc.P;
    mpc.H_f0_ter(mpc.s_col,mpc.s_col) = mpc.H_f0_ter(mpc.s_col,mpc.s_col) + mpc.P2;
end

end

function mpc = genCustomCost(mpc)

s_index = mpc.s_col;
su_index = mpc.su_col;
u_index = mpc.u_col;

% compute gradient of z
for k = 1:mpc.N-1
    if mpc.z_use_s && mpc.z_use_su && mpc.z_use_u
        mpc.grad_z(:,:,k) = [mpc.Cz(:,:,k)'; mpc.Dsuz(:,:,k)'; mpc.Dz(:,:,k)'];
        index_k = [s_index su_index u_index];

    elseif mpc.z_use_s && mpc.z_use_su
        mpc.grad_z(:,:,k) = [mpc.Cz(:,:,k)'; mpc.Dsuz(:,:,k)'];
        index_k = [s_index su_index];

    elseif mpc.z_use_s && mpc.z_use_u
        mpc.grad_z(:,:,k) = [mpc.Cz(:,:,k)'; mpc.Dz(:,:,k)'];
        index_k = [s_index u_index];

    elseif mpc.z_use_su && mpc.z_use_u
        mpc.grad_z(:,:,k) = [mpc.Dsuz(:,:,k)'; mpc.Dz(:,:,k)'];
        index_k = [su_index u_index];

    elseif mpc.z_use_s
        mpc.grad_z(:,:,k) = mpc.Cz(:,:,k)';
        index_k = s_index;

    elseif mpc.z_use_su
        mpc.grad_z(:,:,k) = mpc.Dsuz(:,:,k)';
        index_k = su_index;

    elseif mpc.z_use_u
        mpc.grad_z(:,:,k) = mpc.Dz(:,:,k)';
        index_k = u_index;
    end
end
mpc.custom_cost_index_k = index_k;

if mpc.z_use_k0
    mpc.grad_z_0 = mpc.Dz_0';
end
if mpc.z_use_ter
    mpc.grad_z_ter = mpc.Cz_ter';
end

% gradient/Hessian of quadratic cost on z
if mpc.quad_custom_cost

    % k = 0
    if mpc.z_use_k0
        gradz_Qz_0 = mpc.grad_z_0*mpc.Qz_0;
        H_CustomCost_0 = gradz_Qz_0*mpc.grad_z_0';

        mpc.gradz_Qz_0 = gradz_Qz_0;
        mpc.H_f0_0 = mpc.H_f0_0 + H_CustomCost_0;
        mpc.H_CustomCost_0 = H_CustomCost_0;
    end

    for k = 1:mpc.N-1

        gradz_Qz = mpc.grad_z(:,:,k)*mpc.Qz(:,:,k);
        H_CustomCost = gradz_Qz*mpc.grad_z(:,:,k)';

        mpc.gradz_Qz_k(:,:,k) = gradz_Qz;
        mpc.H_f0_k(index_k,index_k,k) = mpc.H_f0_k(index_k,index_k,k)...
                                        + H_CustomCost;
        mpc.H_CustomCost_k(:,:,k) = H_CustomCost;
    end

    if mpc.z_use_ter
        gradz_Qz_ter = mpc.grad_z_ter*mpc.Qz_ter;
        H_CustomCost_ter = gradz_Qz_ter*mpc.grad_z_ter';

        mpc.gradz_Qz_ter = gradz_Qz_ter;
        mpc.H_f0_ter(mpc.s_col,mpc.s_col) = mpc.H_f0_ter(mpc.s_col,mpc.s_col) +...
            H_CustomCost_ter;
        mpc.H_CustomCost_ter = H_CustomCost_ter;
    end
end

% gradient of linear cost on z
if mpc.lin_custom_cost

    % k = 0
    if mpc.z_use_k0
        mpc.gradz_qz_0 = mpc.grad_z_0*mpc.qz_0;
    end

    for k = 1:mpc.N-1
        mpc.gradz_qz_k(:,k) = mpc.grad_z(:,:,k)*mpc.qz;
    end

    if mpc.z_use_ter
        mpc.gradz_qz_ter = mpc.grad_z_ter*mpc.qz_ter;
    end
end

end

function mpc = genTrackingCost(mpc)

s_index = mpc.s_col;
u_index = mpc.u_col;

if mpc.y_use_k0
    grad_err_0 = -mpc.D_0';
    gradErr_Qe_0 = grad_err_0*mpc.Qe_0;
    H_ErrCost_0 =  gradErr_Qe_0*grad_err_0';

    mpc.grad_err_0 = grad_err_0;
    mpc.gradErr_Qe_0 = gradErr_Qe_0;
    mpc.H_f0_0(:,:) = mpc.H_f0_0(:,:) + H_ErrCost_0;
    mpc.H_ErrCost_0 = H_ErrCost_0;
end

for k = 1:mpc.N-1

    % gradient of err = ref-y
    if mpc.y_use_s && mpc.y_use_u
        grad_err = [-mpc.C(:,:,k)'; -mpc.D(:,:,k)'];
        index_k = [s_index u_index];
    elseif mpc.y_use_s
        grad_err = -mpc.C(:,:,k)';
        index_k = s_index;
    elseif mpc.y_use_u
        grad_err = -mpc.D(:,:,k)';
        index_k = u_index;
    end

    % gradient/hessian of error cost
    gradErr_Qe = grad_err*mpc.Qe(:,:,k);
    H_ErrCost = gradErr_Qe*grad_err';
    
    mpc.grad_err(:,:,k) = grad_err;
    mpc.gradErr_Qe_k(:,:,k) = gradErr_Qe;
    mpc.H_f0_k(index_k,index_k,k) = mpc.H_f0_k(index_k,index_k,k) + ...
                                    + H_ErrCost;
    mpc.H_ErrCost_k(:,:,k) = H_ErrCost;
end
mpc.tracking_cost_index_k = index_k;

if mpc.y_use_ter
    grad_err_ter = -mpc.C_ter';
    gradErr_Qe_ter = grad_err_ter*mpc.Qe_ter;
    H_ErrCost_ter =  gradErr_Qe_ter*grad_err_ter';

    mpc.grad_err_ter = grad_err_ter;
    mpc.gradErr_Qe_ter = gradErr_Qe_ter;
    mpc.H_f0_ter(mpc.s_col,mpc.s_col) = mpc.H_f0_ter(mpc.s_col,mpc.s_col) +...
                                        H_ErrCost_ter;
    mpc.H_ErrCost_ter = H_ErrCost_ter;
end

end

function mpc = genControlCost(mpc)

if mpc.quad_control_cost

    mpc.H_f0_0 = mpc.H_f0_0 + mpc.Ru(:,:,1);

    for k = 1:mpc.N-1
        u_index = mpc.u_col;

        mpc.H_f0_k(u_index,u_index,k) = mpc.H_f0_k(u_index,u_index,k) + ...
                                            + mpc.Ru(:,:,k+1);
    end
end

end

function mpc = genDiffControlCost(mpc)

su_index = mpc.su_col;
u_index = mpc.u_col;

index_k = [su_index u_index];
mpc.du_index_k = index_k;

% k = 0
mpc.H_f0_0 = mpc.H_f0_0 + mpc.Rdu(:,:,1);

for k = 1:mpc.N-1

    Rdu_k = mpc.Rdu(:,:,k+1);
    
    gradRateCtrl_Rdu = [-Rdu_k;Rdu_k];
    H_RateCtrl = [Rdu_k -Rdu_k;-Rdu_k Rdu_k];

    mpc.gradRateCtrl_Rdu_k(:,:,k) = gradRateCtrl_Rdu;
    mpc.H_f0_k(index_k,index_k,k) = mpc.H_f0_k(index_k,index_k,k) + ...
                                        + H_RateCtrl;
    mpc.H_RateCtrl_k(:,:,k) = H_RateCtrl;
end

end

function mpc = genSoftSlacksCost(mpc)

grad_qv_0 = zeros(mpc.nv_k(1),1);
grad_qv_k = zeros(mpc.nv_k(2),mpc.N-1);
grad_qv_ter = zeros(mpc.nv_k(3),1);

if mpc.nv_k(1)
    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit && mpc.y_use_k0
            row = mpc.y_cnstr.min_row_v_0;
            grad_qv_0(row) = mpc.y_cnstr.qv_min_0;
        end
        if mpc.y_cnstr.max_limit && mpc.y_use_k0
            row = mpc.y_cnstr.max_row_v_0;
            grad_qv_0(row) = mpc.y_cnstr.qv_max_0;
        end
    end

    if mpc.has_h_cnstr
        if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_k0
            row = mpc.h_cnstr.min_row_v_0;
            grad_qv_0(row) = mpc.h_cnstr.qv_min_0;
        end
        if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_k0
            row = mpc.h_cnstr.max_row_v_0;
            grad_qv_0(row) = mpc.h_cnstr.qv_max_0;
        end
    end
end

if mpc.nv_k(2)
for k = 1:mpc.N-1
    if mpc.has_s_cnstr
        if mpc.s_cnstr.min_limit
            row = mpc.s_cnstr.min_row_v_k;
            grad_qv_k(row,k) = mpc.s_cnstr.qv_min;
        end
        if mpc.s_cnstr.max_limit
            row = mpc.s_cnstr.max_row_v_k;
            grad_qv_k(row,k) = mpc.s_cnstr.qv_max;
        end
    end

    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit
            row = mpc.y_cnstr.min_row_v_k;
            grad_qv_k(row,k) = mpc.y_cnstr.qv_min;
        end
        if mpc.y_cnstr.max_limit
            row = mpc.y_cnstr.max_row_v_k;
            grad_qv_k(row,k) = mpc.y_cnstr.qv_max;
        end
    end

    if mpc.has_h_cnstr
        if mpc.h_cnstr.min_limit
            row = mpc.h_cnstr.min_row_v_k;
            grad_qv_k(row,k) = mpc.h_cnstr.qv_min;
        end
        if mpc.h_cnstr.max_limit
            row = mpc.h_cnstr.max_row_v_k;
            grad_qv_k(row,k) = mpc.h_cnstr.qv_max;
        end
    end
    
end
end

% k = N

if mpc.ng_k(3)
    
if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit
        row = mpc.s_cnstr.min_ineqRow_ter;
        grad_qv_ter(row) = mpc.s_cnstr.qv_min;
    end
    if mpc.s_cnstr.max_limit
        row = mpc.s_cnstr.max_ineqRow_ter;
        grad_qv_ter(row) = mpc.s_cnstr.qv_max;
    end
end

if mpc.has_y_cnstr
    if mpc.y_cnstr.min_limit && mpc.y_use_ter
        row = mpc.y_cnstr.min_ineqRow_ter;
        grad_qv_ter(row) = mpc.y_cnstr.qv_min_ter;
    end
    if mpc.y_cnstr.max_limit && mpc.y_use_ter
        row = mpc.y_cnstr.max_ineqRow_ter;
        grad_qv_ter(row) = mpc.y_cnstr.qv_max_ter;
    end
end

if mpc.has_h_cnstr
    if mpc.h_cnstr.min_limit && mpc.h_cnstr.use_ter
        row = mpc.h_cnstr.min_ineqRow_ter;
        grad_qv_ter(row) = mpc.h_cnstr.qv_min_ter;
    end
    if mpc.h_cnstr.max_limit && mpc.h_cnstr.use_ter
        row = mpc.h_cnstr.max_ineqRow_ter;
        grad_qv_ter(row) = mpc.h_cnstr.qv_max_ter;
    end
end

end

mpc.grad_qv_0 = grad_qv_0;
mpc.grad_qv_k = grad_qv_k;
mpc.grad_qv_ter = grad_qv_ter;

end