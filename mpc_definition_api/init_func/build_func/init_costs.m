function mpc = init_costs(mpc)

mpc.hessCost = zeros(mpc.n);

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
    s_index = mpc.s_index_k(:,mpc.N+1);

    mpc.P2 = -2*mpc.P;
    mpc.hessCost(s_index,s_index) = ...
            mpc.hessCost(s_index,s_index) + 2*mpc.P;
end

end

function mpc = genCustomCost(mpc)

if mpc.quad_custom_cost

    % k = 0
    if mpc.z_use_u
        gradPerfQz = mpc.Dz'*mpc.Qz;
        hessPerfCost = gradPerfQz*mpc.Dz;

        u_index = mpc.u_index_k(:,1);

        mpc.gradPerfQz_0 = gradPerfQz;
        mpc.hessCost(u_index,u_index) = ...
            mpc.hessCost(u_index,u_index) + hessPerfCost;
    end

    for k = 2:mpc.N

        index_k = [];
        s_index = mpc.s_index_k(:,k);
        su_index = mpc.su_index_k(:,k);
        u_index = mpc.u_index_k(:,k);

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

        gradPerfQz = grad_z*mpc.Qz;
        hessPerfCost = gradPerfQz*grad_z';

        mpc.gradPerfQz(:,:,k) = gradPerfQz;
        mpc.hessCost(index_k,index_k) = ...
            mpc.hessCost(index_k,index_k) + hessPerfCost;
    end
end

if mpc.lin_custom_cost

    mpc.gradPerfqz = zeros(mpc.n,1);

    % k = 0
    if mpc.z_use_u
        u_index = mpc.u_index_k(:,1);
        mpc.gradPerfqz(u_index) = mpc.Dz'*mpc.qz;
    end

    for k = 2:mpc.N

        index_k = [];
        s_index = mpc.s_index_k(:,k);
        su_index = mpc.su_index_k(:,k);
        u_index = mpc.u_index_k(:,k);

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

        mpc.gradPerfqz(index_k) = grad_z*mpc.qz;
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

for k = 2:mpc.N

    index_k = [];
    s_index = mpc.s_index_k(:,k);
    u_index = mpc.u_index_k(:,k);

    grad_err = [];

    if use_s
        index_k = [index_k;s_index];
        grad_err =  [grad_err;-mpc.C'];
    end
    if use_u
        index_k = [index_k;u_index];
        grad_err =  [grad_err;-mpc.D'];
    end

    gradErrQe = grad_err*mpc.Qe;
    hessErrCost = gradErrQe*grad_err';
    
    mpc.gradErrQe(:,:,k) = gradErrQe;
    mpc.hessCost(index_k,index_k) = ...
        mpc.hessCost(index_k,index_k) + hessErrCost;
end
end

function mpc = genControlCost(mpc)

if mpc.quad_control_cost

    for k = 1:mpc.N
        u_index = mpc.u_index_k(:,k);

        mpc.gradCtlrRu(:,:,k) = mpc.Ru;
        mpc.hessCost(u_index,u_index) = ...
            mpc.hessCost(u_index,u_index) + mpc.Ru;
    end
end

if mpc.lin_control_cost
    mpc.gradCtlrru = zeros(mpc.n,1);
    for k = 1:mpc.N
        u_index = mpc.u_index_k(:,k);
        mpc.gradCtlrru(u_index) = mpc.ru;
    end
end

end

function mpc = genDiffControlCost(mpc)
% k = 0
u_index = mpc.u_index_k(:,1);

mpc.gradDiffCtlrR_0 = mpc.Rdu;
mpc.hessCost(u_index,u_index) = ...
    mpc.hessCost(u_index,u_index) + mpc.Rdu;

for k = 2:mpc.N
    su_index = mpc.su_index_k(:,k);
    u_index = mpc.u_index_k(:,k);

    index_k = [su_index;u_index];

    gradDUCost = [-eye(mpc.nu);eye(mpc.nu)]*mpc.Rdu;
    hessDUCost = gradDUCost*[-eye(mpc.nu) eye(mpc.nu)];

    mpc.gradDiffCtlrR(:,:,k) = gradDUCost;
    mpc.hessCost(index_k,index_k) = ...
        mpc.hessCost(index_k,index_k) + hessDUCost;
end

end

function mpc = genSoftSlacksCost(mpc)
gradSlackqv = zeros(mpc.n,1);

for k = 1:mpc.N+1

    if mpc.has_s_cnstr
        if mpc.s_cnstr.min_limit && any(mpc.s_cnstr.v_min_index_k(:,k))
            gradSlackqv(mpc.s_cnstr.v_min_index_k(:,k)) = mpc.s_cnstr.qv_min;
        end
        if mpc.s_cnstr.max_limit && any(mpc.s_cnstr.v_max_index_k(:,k))
            gradSlackqv(mpc.s_cnstr.v_max_index_k(:,k)) = mpc.s_cnstr.qv_max;
        end
    end
    if mpc.has_y_cnstr
        if mpc.y_cnstr.min_limit && any(mpc.y_cnstr.v_min_index_k(:,k))
            gradSlackqv(mpc.y_cnstr.v_min_index_k(:,k)) = mpc.y_cnstr.qv_min;
        end
        if mpc.y_cnstr.max_limit && any(mpc.y_cnstr.v_max_index_k(:,k))
            gradSlackqv(mpc.y_cnstr.v_max_index_k(:,k)) = mpc.y_cnstr.qv_max;
        end
    end
    if mpc.has_h_cnstr
        if mpc.h_cnstr.min_limit && any(mpc.h_cnstr.v_min_index_k(:,k))
            gradSlackqv(mpc.h_cnstr.v_min_index_k(:,k)) = mpc.h_cnstr.qv_min;
        end
        if mpc.h_cnstr.max_limit && any(mpc.h_cnstr.v_max_index_k(:,k))
            gradSlackqv(mpc.h_cnstr.v_max_index_k(:,k)) = mpc.h_cnstr.qv_max;
        end
    end
end
mpc.gradSlackqv = gradSlackqv;

end