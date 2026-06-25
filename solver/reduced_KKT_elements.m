function mpc = reduced_KKT_elements(mpc)

%% ri hat
% k = 0
mpc.ri_hat_0 = mpc.ri_0+mpc.g_0;
if any(mpc.v_index_0)
    mpc.ri_hat_0(mpc.vi_0) = mpc.ri_hat_0(mpc.vi_0) - mpc.v_0 +...
        (mpc.v_0.^2).*mpc.grad_qv_0*mpc.t;
end

for k = 1:mpc.N-1
    mpc.ri_hat_k(:,k) = mpc.ri_k(:,k) + mpc.g_k(:,k);
    if mpc.nv_k(2)
        mpc.ri_hat_k(mpc.vi_k,k) = mpc.ri_hat_k(mpc.vi_k,k) - ...
            mpc.v_k(:,k) + (mpc.v_k(:,k).^2).*mpc.grad_qv_k(:,k)*mpc.t;
    end
end

if mpc.has_s_cnstr
    mpc.ri_hat_ter = mpc.ri_ter + mpc.g_ter - mpc.v_ter...
          + (mpc.v_ter.^2).*mpc.grad_qv_ter*mpc.t;
end
%% S

mpc.S_0(:) = mpc.g_0.^2;
if mpc.nv_k(1)
    mpc.S_0(:) = mpc.S_0 + mpc.v_0.^2;
end
mpc.iS_0(:) = 1./mpc.S_0;

mpc.S_k(:,:) = mpc.g_k.^2;
if mpc.nv_k(2)
    mpc.S_k(mpc.vi_k,:) = mpc.S_k(mpc.vi_k,:) + mpc.v_k.^2;
end
mpc.iS_k(:,:) = 1./mpc.S_k;

if mpc.ng_k(mpc.N+1)
    mpc.S_ter(:) = mpc.g_ter.^2+mpc.v_ter.^2;
    mpc.iS_ter(:) = 1./mpc.S_ter;
end

%% rx hat

rx_ineq_0 = mpc.Ai_0'*(mpc.ri_hat_0./mpc.iS_0);
for k = 1:mpc.N-1
    rx_ineq_k(:,k) = mpc.Ai_k(:,:,k)'*(mpc.ri_hat_k(:,k)./mpc.iS_k(:,k));
end
if mpc.ng_k(mpc.N*1)
    rx_ineq_ter = mpc.Ai_ter'*(mpc.ri_hat_ter./mpc.iS_ter);
end

mpc.ru_hat_k(1) = mpc.ru_k(1) +  rx_ineq_0;

mpc.rse_hat_k(1:mpc.nse,:) = mpc.rse_k + rx_ineq_k(1:mpc.nse,:);
mpc.ru_hat_k(2:mpc.N) = mpc.ru_k(2:mpc.N) +  rx_ineq_k(mpc.nse+1:mpc.nse+mpc.nu,:);

mpc.rs_hat_ter(:) = mpc.rs_ter + rx_ineq_ter;

%% H

mpc.R_k(:,:,1) = mpc.t*mpc.H_f0_0 + mpc.Ai_0'*diag(mpc.iS_0)*mpc.Ai_0;
for k = 1:mpc.N-1
    mpc.H_k(:,:,k) = mpc.t*mpc.H_f0_k(:,:,k) + ...
                        mpc.Ai_k(:,:,k)'*diag(mpc.iS_k(:,k))*mpc.Ai_k(:,:,k);

    mpc.Q_k(:,:,k) = mpc.H_k(1:mpc.nse,1:mpc.nse,k);
    mpc.R_k(:,:,k+1) = mpc.H_k(mpc.nse+1:mpc.nvar_k,mpc.nse+1:mpc.nvar_k);
    mpc.Y_k(:,:,k) = mpc.H_k(mpc.nse+1:mpc.nvar_k,1:mpc.nse,k);
end

mpc.Q_ter = mpc.t*mpc.H_f0_ter;
if mpc.ng_k(mpc.N*1)
    mpc.Q_ter(:,:) = mpc.Q_ter + mpc.Ai_ter'*diag(mpc.iS_ter)*mpc.Ai_ter;
end

end