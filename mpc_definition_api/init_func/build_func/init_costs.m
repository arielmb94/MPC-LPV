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

if mpc.diffcontrol_cost
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

index_k = [];
s_index = mpc.s_col;
su_index = mpc.su_col;
u_index = mpc.u_col;

grad_z = [];

if mpc.z_use_s
    index_k = [index_k;s_index];
    grad_z =  [grad_z;mpc.Cz'];
end
if mpc.z_use_su
    index_k = [index_k;su_index];
    grad_z =  [grad_z;mpc.Dsuz'];
end
if mpc.z_use_u
    index_k = [index_k;u_index];
    grad_z =  [grad_z;mpc.Dz'];
end

mpc.grad_z = grad_z;
mpc.custom_cost_index_k = index_k;

if mpc.quad_custom_cost

    % k = 0
    if mpc.z_use_k0
        grad_z_0 = mpc.Dz_0';
        gradPerfQz_0 = grad_z_0*mpc.Qz_0;
        hessPerfCost_0 = gradPerfQz_0*grad_z_0';

        mpc.grad_z_0 = grad_z_0;
        mpc.gradPerfQz_0 = gradPerfQz_0;
        mpc.H_f0_0 = mpc.H_f0_0 + hessPerfCost_0;
    end

    for k = 1:mpc.N-1

        gradPerfQz = grad_z*mpc.Qz;
        hessPerfCost = gradPerfQz*grad_z';

        mpc.gradPerfQz_k(:,:,k) = gradPerfQz;
        mpc.H_f0_k(index_k,index_k,k) = mpc.H_f0_k(index_k,index_k,k) + ...
            + hessPerfCost;
    end

    if mpc.z_use_ter
        grad_z_ter = mpc.Cz_ter';
        gradPerfQz_ter = grad_z_ter*mpc.Qz_ter;
        hessPerfCost_ter = gradPerfQz_ter*grad_z_ter';

        mpc.grad_z_ter = grad_z_ter;
        mpc.gradPerfQz_ter = gradPerfQz_ter;
        mpc.H_f0_ter(mpc.s_col,mpc.s_col) = mpc.H_f0_ter(mpc.s_col,mpc.s_col) +...
            hessPerfCost_ter;
    end
end

if mpc.lin_custom_cost

    % k = 0
    if mpc.z_use_k0
        mpc.gradPerfqz_0 = mpc.Dz'*mpc.qz_0;
    end

    for k = 1:mpc.N-1
        mpc.gradPerfqz_k(:,k) = grad_z*mpc.qz;
    end

    if mpc.z_use_ter
        mpc.gradPerfqz_ter = mpc.Cz_ter'*mpc.qz_ter;
    end
end

end

function mpc = genTrackingCost(mpc)

index_k = [];
s_index = mpc.s_col;
u_index = mpc.u_col;

grad_err = [];

if mpc.y_use_s
    index_k = [index_k;s_index];
    grad_err =  [grad_err;-mpc.C'];
end
if mpc.y_use_u
    index_k = [index_k;u_index];
    grad_err =  [grad_err;-mpc.D'];
end

mpc.tracking_cost_index_k = index_k;
mpc.grad_err = grad_err;

if mpc.y_use_k0
    grad_err_0 = -mpc.D_0';
    gradErrQe_0 = grad_err_0*mpc.Qe_0;
    hessErrCost_0 =  gradErrQe_0*grad_err_0';

    mpc.grad_err_0 = grad_err_0;
    mpc.gradErrQe_0 = gradErrQe_0;
    mpc.H_f0_0(:,:) = mpc.H_f0_0(:,:) + hessErrCost_0;
end

for k = 1:mpc.N-1

    gradErrQe = grad_err*mpc.Qe;
    hessErrCost = gradErrQe*grad_err';
    
    mpc.gradErrQe_k(:,:,k) = gradErrQe;
    mpc.H_f0_k(index_k,index_k,k) = mpc.H_f0_k(index_k,index_k,k) + ...
                                    + hessErrCost;
end

if mpc.y_use_ter
    grad_err_ter = -mpc.C_ter';
    gradErrQe_ter = grad_err_ter*mpc.Qe_ter;
    hessErrCost_ter =  gradErrQe_ter*grad_err_ter';

    mpc.grad_err_ter = grad_err_ter;
    mpc.gradErrQe_ter = gradErrQe_ter;
    mpc.H_f0_ter(mpc.s_col,mpc.s_col) = mpc.H_f0_ter(mpc.s_col,mpc.s_col) +...
                                        hessErrCost_ter;
end

end

function mpc = genControlCost(mpc)

if mpc.quad_control_cost

    mpc.gradCtlrRu_0 = mpc.Ru;
    mpc.H_f0_0 = mpc.H_f0_0 + mpc.Ru;

    for k = 1:mpc.N-1
        u_index = mpc.u_col;

        mpc.gradCtlrRu_k(:,:,k) = mpc.Ru;
        mpc.H_f0_k(u_index,u_index,k) = mpc.H_f0_k(u_index,u_index,k) + ...
                                            + mpc.Ru;
    end
end

if mpc.lin_control_cost
    mpc.gradCtlrru_0 = mpc.ru;
    for k = 1:mpc.N-1
        mpc.gradCtlrru_k(:,k) = mpc.ru;
    end
end

end

function mpc = genDiffControlCost(mpc)
% k = 0
mpc.gradDiffCtlrR_0 = mpc.Rdu;
mpc.H_f0_0 = mpc.H_f0_0 + mpc.Rdu;


su_index = mpc.su_col;
u_index = mpc.u_col;

index_k = [su_index;u_index];
mpc.du_index_k = index_k;

grad_du = [-eye(mpc.nu);eye(mpc.nu)];
mpc.grad_du = grad_du;

for k = 1:mpc.N-1
    gradDUCost = grad_du*mpc.Rdu;
    hessDUCost = gradDUCost*grad_du';

    mpc.gradDiffCtlrR_k(:,:,k) = gradDUCost;
    mpc.H_f0_k(index_k,index_k,k) = mpc.H_f0_k(index_k,index_k,k) + ...
                                        + hessDUCost;
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