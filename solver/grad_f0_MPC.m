function grad_J = grad_f0_MPC(mpc)

grad_J = zeros(mpc.n,1);

if any(mpc.nv_k)
    grad_J(mpc.v_index) = mpc.gradSlackqv(mpc.v_index);
end

% k=0
if mpc.quad_control_cost
    index = mpc.u_index_k(:,1);
    grad_J(index) = grad_J(index) + mpc.gradCtlrRu(:,:,1)*mpc.u(:,1);
end
if mpc.lin_control_cost
    index = mpc.u_index_k(:,1);
    grad_J(index) = grad_J(index) + mpc.gradCtlrru(index);
end

if mpc.diffcontrol_cost
    index = mpc.u_index_k(:,1);
    grad_J(index) = grad_J(index) + mpc.gradDiffCtlrR_0*mpc.du(:,1);
end

if mpc.quad_custom_cost
    if mpc.z_use_u
        index = mpc.u_index_k(:,1);
        grad_J(index) = grad_J(index) + mpc.gradPerfQz_0*mpc.z(:,1);
    end
end
if mpc.lin_custom_cost
    if mpc.z_use_u
        index = mpc.u_index_k(:,1);
        grad_J(index) = grad_J(index) + mpc.gradPerfqz(index);
    end
end

for k = 2:mpc.N

    if mpc.tracking_cost
        index = [];
        if mpc.y_use_s
            index = [index;mpc.s_index_k(:,k)];
        end
        if mpc.y_use_u
            index = [index;mpc.u_index_k(:,k)];
        end
        grad_J(index) = grad_J(index) + mpc.gradErrQe(:,:,k)*mpc.err(:,k-1);
    end

    if mpc.quad_control_cost
        index = mpc.u_index_k(:,k);
        grad_J(index) = grad_J(index) + mpc.gradCtlrRu(:,:,k)*mpc.u(:,k);
    end
    if mpc.lin_control_cost
        index = mpc.u_index_k(:,k);
        grad_J(index) = grad_J(index) + mpc.gradCtlrru(index);
    end

    if mpc.diffcontrol_cost
        index = [mpc.su_index_k(:,k); mpc.u_index_k(:,k)];
        grad_J(index) = grad_J(index) + mpc.gradDiffCtlrR(:,:,k)*mpc.du(:,k);
    end

    if mpc.quad_custom_cost || mpc.lin_custom_cost
        index = [];
        if mpc.z_use_s
            index = [index;mpc.s_index_k(:,k)];
        end
        if mpc.z_use_su
            index = [index;mpc.su_index_k(:,k)];
        end
        if mpc.z_use_u
            index = [index;mpc.u_index_k(:,k)];
        end

        if mpc.quad_custom_cost
            grad_J(index) = grad_J(index) + mpc.gradPerfQz(:,:,k)*mpc.z(:,k);
        end
        if mpc.lin_custom_cost
            grad_J(index) = grad_J(index) + mpc.gradPerfqz(index);
        end
    end
end

if mpc.ter_ingredients
    index = mpc.s_index_k(:,mpc.N+1);
    grad_J(index) = grad_J(index) - mpc.P2*(mpc.xN_ref-mpc.s_ter);
end
end