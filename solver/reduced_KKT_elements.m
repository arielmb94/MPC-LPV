function mpc = reduced_KKT_elements(mpc)
%% ri_hat

mpc.ri_hat_0 = mpc.ri_0+mpc.g_0;
if mpc.nv_k(1)
    mpc.rv_0 = mpc.t*mpc.grad_qv_0.*(mpc.v_0.^2)-mpc.v_0;
    mpc.ri_hat_0 = mpc.ri_hat_0(mpc.v_rows_0) + mpc.rv_0;
end

mpc.ri_hat_k = mpc.ri_k+mpc.g_k;
if mpc.nv_k(2)
    mpc.rv_k = mpc.t*mpc.grad_qv_k.*(mpc.v_k.^2)-mpc.v_k;
    mpc.ri_hat_k(mpc.v_rows_k,:) = mpc.ri_hat_k(mpc.v_rows_k,:) + mpc.rv_k;
end

if mpc.ng_k(3)
    mpc.rv_ter = mpc.t*mpc.grad_qv_ter.*(mpc.v_ter.^2)-mpc.v_ter;
    mpc.ri_hat_ter = mpc.ri_ter+mpc.g_ter + mpc.rv_ter;
end

%% S

mpc.iS_0 = mpc.g_0.^2;
if mpc.nv_k(1), mpc.iS_0(mpc.v_rows_0) = mpc.iS_0(mpc.v_rows_0) + mpc.v_0.^2; end
mpc.iS_0  = 1./mpc.iS_0;

mpc.iS_k = mpc.g_k.^2;
if mpc.nv_k(2), mpc.iS_k(mpc.v_rows_k,:) = mpc.iS_k(mpc.v_rows_k,:) + mpc.v_k.^2; end
mpc.iS_k  = 1./mpc.iS_k;

mpc.iS_ter = mpc.g_ter.^2;
if mpc.nv_k(3), mpc.iS_ter = mpc.iS_ter + mpc.v_ter.^2; end
mpc.iS_ter  = 1./mpc.iS_ter;

%% rx_hat

mpc.iS_ri_hat_0 = mpc.iS_0.*mpc.ri_hat_0;
mpc.ru_hat_0 = mpc.ru_0 + mpc.Ai_0'*mpc.iS_ri_hat_0;


mpc.iS_ri_hat_k = mpc.iS_k.*mpc.ri_hat_k;
for k = 1:mpc.N-1

    mpc.rse_hat_k(:,k) = mpc.rse_k(:,k) + mpc.Ai_k(:,mpc.se_col,k)'*mpc.iS_ri_hat_k(:,k);
    mpc.ru_hat_k(:,k) = mpc.ru_k(:,k) + mpc.Ai_k(:,mpc.u_col,k)'*mpc.iS_ri_hat_k(:,k);
end

mpc.iS_ri_hat_ter = mpc.iS_ter.*mpc.ri_hat_ter;
mpc.rse_hat_ter = mpc.rse_ter + mpc.Ai_ter'*mpc.iS_ri_hat_ter;

%% H = Hess(f0) + Ai'*(S^-1)*Ai

mpc.iS_Ai_0 = mpc.iS_0.*mpc.Ai_0;
mpc.R_0(:,:) = mpc.t*mpc.H_f0_0 + mpc.Ai_0'*mpc.iS_Ai_0;

for k = 1:mpc.N-1

    mpc.iS_Ai_k(:,:,k) = mpc.iS_k(:,k).*mpc.Ai_k(:,:,k);

    mpc.H_k(:,:,k) = mpc.t*mpc.H_f0_k(:,:,k) + ...
                        mpc.Ai_k(:,:,k)'*mpc.iS_Ai_k(:,:,k);

    mpc.Q_k(:,:,k) = mpc.H_k(mpc.se_col,mpc.se_col,k);
    mpc.R_k(:,:,k) = mpc.H_k(mpc.u_col,mpc.u_col,k);
    mpc.Y_k(:,:,k) = mpc.H_k(mpc.u_col,mpc.se_col,k);
end

mpc.Q_ter(:,:) = mpc.t*mpc.H_f0_ter;
if mpc.ng_k(3)
    mpc.iS_Ai_ter = mpc.iS_ter.*mpc.Ai_ter;
    mpc.Q_ter(:,:) = mpc.Q_ter(:,:) + mpc.Ai_ter'*mpc.iS_Ai_ter;
end

end
