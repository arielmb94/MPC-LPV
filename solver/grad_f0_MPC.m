function mpc = grad_f0_MPC(mpc)

mpc.grad_f0_0(:) = 0;
mpc.grad_f0_k(:,:) = 0;
mpc.grad_f0_ter(:) = 0;

% k=0
if mpc.quad_control_cost
    mpc.grad_f0_0 = mpc.grad_f0_0 + mpc.gradCtlrRu(:,:,1)*mpc.u(:,1);
end
if mpc.lin_control_cost
    mpc.grad_f0_0 = mpc.grad_f0_0 + mpc.gradCtlrru(:,1);
end

if mpc.diffcontrol_cost
    mpc.grad_f0_0 = mpc.grad_f0_0 + mpc.gradDiffCtlrR_0*mpc.du(:,1);
end

if mpc.quad_custom_cost
    if mpc.z_use_u
        mpc.grad_f0_0 = mpc.grad_f0_0 + mpc.gradPerfQz_0*mpc.z(:,1);
    end
end
if mpc.lin_custom_cost
    if mpc.z_use_u
        mpc.grad_f0_0 = mpc.grad_f0_0 + mpc.mpc.gradPerfqz_0;
    end
end

for k = 1:mpc.N-1

    if mpc.tracking_cost
        index = mpc.tracking_cost_index_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradErrQe(:,:,k)*mpc.err(:,k);
    end

    if mpc.quad_control_cost
        index = mpc.nse+1:mpc.nvar_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradCtlrRu(:,:,k+1)*mpc.u(:,k+1);
    end
    if mpc.lin_control_cost
        index = mpc.nse+1:mpc.nvar_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradCtlrru(:,k+1);
    end

    if mpc.diffcontrol_cost
        index = mpc.du_index_k;
        mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradDiffCtlrR(:,:,k)*mpc.du(:,k+1);
    end

    if mpc.quad_custom_cost || mpc.lin_custom_cost
        index = mpc.custom_cost_index_k;

        if mpc.quad_custom_cost
            mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradPerfQz(:,:,k)*mpc.z(:,k+1);
        end
        if mpc.lin_custom_cost
            mpc.grad_f0_k(index,k) = mpc.grad_f0_k(index,k) + mpc.gradPerfqz_k(:,k);
        end
    end
end

if mpc.ter_ingredients
    mpc.grad_f0_ter(1:mpc.nx) = mpc.grad_f0_ter(1:mpc.nx) - mpc.P2*(mpc.xN_ref-mpc.s_ter);
end

u_index = mpc.nse+1:mpc.nvar_k;
mpc.ru_k(:,:) = mpc.t*[mpc.grad_f0_0 mpc.grad_f0_k(u_index,:)];
mpc.rse_k(:,1:mpc.N-1) = mpc.t*mpc.grad_f0_k(1:mpc.nse,:);
mpc.rse_k(:,mpc.N) = mpc.t*mpc.grad_f0_ter;

end