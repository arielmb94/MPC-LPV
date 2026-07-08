function mpc = init_costs(mpc)

mpc.nvar_k = mpc.nse+mpc.nu;
mpc.H_f0_0 = zeros(mpc.nu);
mpc.H_f0_k = zeros(mpc.nvar_k,mpc.nvar_k,mpc.N-1);
mpc.H_f0_ter = zeros(mpc.nse);

mpc.H_k = zeros(mpc.nvar_k,mpc.nvar_k,mpc.N-1);
mpc.Q_k = zeros(mpc.nse,mpc.nse,mpc.N);
mpc.R_k = zeros(mpc.nu,mpc.nu,mpc.N);
mpc.Y_k = zeros(mpc.nu,mpc.nse,mpc.N-1);

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
    mpc.H_f0_ter(1:mpc.nx,1:mpc.nx) = mpc.H_f0_ter(1:mpc.nx,1:mpc.nx) + mpc.P2;
end

end

function mpc = genCustomCost(mpc)

index_k = [];
s_index = 1:mpc.nx;
su_index = mpc.nx+1:mpc.nse;
u_index = mpc.nse+1:mpc.nvar_k;

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
    if mpc.z_use_u
        mpc.grad_z_0 = mpc.Dz';
        gradPerfQz = mpc.Dz'*mpc.Qz;
        hessPerfCost = gradPerfQz*mpc.Dz;

        mpc.gradPerfQz_0 = gradPerfQz;
        mpc.H_f0_0 = mpc.H_f0_0 + hessPerfCost;
    end

    for k = 1:mpc.N-1

        gradPerfQz = grad_z*mpc.Qz;
        hessPerfCost = gradPerfQz*grad_z';

        mpc.gradPerfQz(:,:,k) = gradPerfQz;
        mpc.H_f0_k(index_k,index_k,k) = mpc.H_f0_k(index_k,index_k,k) + ...
                                        + hessPerfCost;
    end
end

if mpc.lin_custom_cost

    mpc.gradPerfqz_0 = zeros(mpc.nu,1);
    mpc.gradPerfqz_k = zeros(length(index_k),mpc.N-1);

    % k = 0
    if mpc.z_use_u
        mpc.gradPerfqz_0 = mpc.Dz'*mpc.qz;
    end

    for k = 1:mpc.N-1
        mpc.gradPerfqz_k(:,k) = grad_z*mpc.qz;
    end
end


end

function mpc = genTrackingCost(mpc)

if ~isempty(mpc.C) && max(any(mpc.C))
    use_s = 1;
else
    use_s = 0;
end
if ~isempty(mpc.D) && max(any(mpc.D))
    use_u = 1;
else
    use_u = 0;
end

index_k = [];
s_index = 1:mpc.nx;
u_index = mpc.nse+1:mpc.nvar_k;

grad_err = [];

if use_s
    index_k = [index_k;s_index];
    grad_err =  [grad_err;-mpc.C'];
end
if use_u
    index_k = [index_k;u_index];
    grad_err =  [grad_err;-mpc.D'];
end

mpc.tracking_cost_index_k = index_k;
mpc.grad_err = grad_err;

for k = 1:mpc.N-1

    gradErrQe = grad_err*mpc.Qe;
    hessErrCost = gradErrQe*grad_err';
    
    mpc.gradErrQe(:,:,k) = gradErrQe;
    mpc.H_f0_k(index_k,index_k,k) = mpc.H_f0_k(index_k,index_k,k) + ...
                                    + hessErrCost;
end
end

function mpc = genControlCost(mpc)

if mpc.quad_control_cost

    mpc.gradCtlrRu(:,:,1) = mpc.Ru;
    mpc.H_f0_0 = mpc.H_f0_0 + mpc.Ru;

    for k = 2:mpc.N
        u_index = mpc.nse+1:mpc.nvar_k;

        mpc.gradCtlrRu(:,:,k) = mpc.Ru;
        mpc.H_f0_k(u_index,u_index,k-1) = mpc.H_f0_k(u_index,u_index,k-1) + ...
                                            + mpc.Ru;
    end
end

if mpc.lin_control_cost
    mpc.gradCtlrru = zeros(mpc.nu,mpc.N);
    for k = 1:mpc.N
        mpc.gradCtlrru(:,k) = mpc.ru;
    end
end

end

function mpc = genDiffControlCost(mpc)
% k = 0
mpc.gradDiffCtlrR_0 = mpc.Rdu;
mpc.H_f0_0 = mpc.H_f0_0 + mpc.Rdu;


su_index = mpc.nx+1:mpc.nse;
u_index = mpc.nse+1:mpc.nvar_k;

index_k = [su_index;u_index];
mpc.du_index_k = index_k;

grad_du = [-eye(mpc.nu);eye(mpc.nu)];
mpc.grad_du = grad_du;

for k = 1:mpc.N-1
    gradDUCost = grad_du*mpc.Rdu;
    hessDUCost = gradDUCost*grad_du';

    mpc.gradDiffCtlrR(:,:,k) = gradDUCost;
    mpc.H_f0_k(index_k,index_k,k) = mpc.H_f0_k(index_k,index_k,k) + ...
                                        + hessDUCost;
end

end

function mpc = genSoftSlacksCost(mpc)

grad_qv_0 = zeros(mpc.nv_k(1),1);
grad_qv_k = zeros(mpc.nv_k(2),1);
grad_qv_ter = zeros(mpc.nv_k(mpc.N+1),1);

if ~isempty(grad_qv_0)
    if mpc.h_cnstr.min_limit
        row = mpc.h_cnstr.min_row_v_0;
        grad_qv_0(row) = mpc.h_cnstr.qv_min;
    end
    if mpc.h_cnstr.max_limit
        row = mpc.h_cnstr.max_row_v_0;
        grad_qv_0(row) = mpc.h_cnstr.qv_max;
    end
end

if ~isempty(grad_qv_k)
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

if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit
        row = mpc.s_cnstr.min_row_ter;
        grad_qv_ter(row) = mpc.s_cnstr.qv_min;
    end
    if mpc.s_cnstr.max_limit
        row = mpc.s_cnstr.max_row_ter;
        grad_qv_ter(row) = mpc.s_cnstr.qv_max;
    end
end

mpc.grad_qv_0 = grad_qv_0;
mpc.grad_qv_k = grad_qv_k;
mpc.grad_qv_ter = grad_qv_ter;

end