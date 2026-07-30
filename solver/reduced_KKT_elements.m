function mpc = reduced_KKT_elements(mpc)
%% rg, rv

mpc.rg_0(:) = -1./(mpc.g_0-mpc.slack_epsilon);
mpc.rg_k(:,:) = -1./(mpc.g_k-mpc.slack_epsilon);
mpc.rg_ter(:) = -1./(mpc.g_ter-mpc.slack_epsilon);

mpc.g2_0 = (mpc.g_0-mpc.slack_epsilon).^2;
mpc.g2_k(:,:) = (mpc.g_k-mpc.slack_epsilon).^2;
mpc.g2_ter(:) = (mpc.g_ter-mpc.slack_epsilon).^2;

if mpc.nv_k(1)
    mpc.rv_0(:) = mpc.t*(mpc.grad_qv_0/mpc.Qv_fctr).*(mpc.v_0-mpc.slack_epsilon) + ...
                  mpc.t*mpc.grad_qv_0 - ...
                  1./(mpc.v_0-mpc.slack_epsilon);
    mpc.v2_0(:) = 1./(mpc.t*mpc.grad_qv_0/mpc.Qv_fctr + 1./((mpc.v_0-mpc.slack_epsilon).^2));
end

if mpc.nv_k(2)
    mpc.rv_k(:,:) = mpc.t*(mpc.grad_qv_k/mpc.Qv_fctr).*(mpc.v_k-mpc.slack_epsilon) + ...
                    mpc.t*mpc.grad_qv_k - ...
                    1./(mpc.v_k-mpc.slack_epsilon);
    mpc.v2_k(:,:) = 1./(mpc.t*mpc.grad_qv_k/mpc.Qv_fctr+1./((mpc.v_k-mpc.slack_epsilon).^2));
end

if mpc.nv_k(3)
    mpc.rv_ter(:) = mpc.t*(mpc.grad_qv_ter/mpc.Qv_fctr).*(mpc.v_ter-mpc.slack_epsilon) + ...
                    mpc.t*mpc.grad_qv_ter - ...
                    1./(mpc.v_ter-mpc.slack_epsilon);
    mpc.v2_ter(:) = 1./(mpc.t*mpc.grad_qv_ter/mpc.Qv_fctr + 1./((mpc.v_ter-mpc.slack_epsilon).^2));
end

%% ri hat = ri - g^2*rg + v^2*rv

% k = 0
mpc.ri_hat_0 = mpc.ri_0 - mpc.g2_0.*mpc.rg_0;
if mpc.nv_k(1)
    mpc.ri_hat_0(mpc.v_rows_0) = mpc.ri_hat_0(mpc.v_rows_0) + mpc.v2_0.*mpc.rv_0;
end

% k = 1,...,N-1
mpc.ri_hat_k(:,:) = mpc.ri_k - mpc.g2_k.*mpc.rg_k;
if mpc.nv_k(2)
    mpc.ri_hat_k(mpc.v_rows_k,:) = mpc.ri_hat_k(mpc.v_rows_k,:) +  mpc.v2_k.*mpc.rv_k;
end

% k = N
if mpc.nv_k(3)
    mpc.ri_hat_ter = mpc.ri_ter - mpc.g2_ter.*mpc.rg_ter + mpc.v2_ter.*mpc.rv_ter;
end
%% S = g^2 + v^2

mpc.S_0(:) = mpc.g2_0;
if mpc.nv_k(1)
    mpc.S_0(mpc.v_rows_0) = mpc.S_0(mpc.v_rows_0) + mpc.v2_0;
end
mpc.iS_0(:) = 1./mpc.S_0;

mpc.S_k(:,:) = mpc.g2_k;
if mpc.nv_k(2)
    mpc.S_k(mpc.v_rows_k,:) = mpc.S_k(mpc.v_rows_k,:) + mpc.v2_k;
end
mpc.iS_k(:,:) = 1./mpc.S_k;

if mpc.ng_k(3)
    mpc.S_ter(:) = mpc.g2_ter + mpc.v2_ter;
    mpc.iS_ter(:) = 1./mpc.S_ter;
end

%% rx hat = rx + Ai'*(S^-1)*ri_hat

mpc.rx_ineq_0(:) = mpc.Ai_0'*(mpc.iS_0.*mpc.ri_hat_0);
for k = 1:mpc.N-1
    mpc.rx_ineq_k(:,k) = mpc.Ai_k(:,:,k)'*(mpc.iS_k(:,k).*mpc.ri_hat_k(:,k));
end
if mpc.ng_k(3)
    mpc.rx_ineq_ter(:) = mpc.Ai_ter'*(mpc.iS_ter.*mpc.ri_hat_ter);
end

mpc.ru_hat_0(:) = mpc.ru_0 +  mpc.rx_ineq_0;

mpc.rse_hat_k(:,:) = mpc.rse_k + mpc.rx_ineq_k(mpc.se_col,:);
mpc.ru_hat_k(:,:) = mpc.ru_k +  mpc.rx_ineq_k(mpc.u_col,:);

mpc.rse_hat_ter(:) = mpc.rse_ter + mpc.rx_ineq_ter;

%% H = Hess(f0) + Ai'*(S^-1)*Ai

mpc.R_0(:,:) = mpc.t*mpc.H_f0_0 + mpc.Ai_0'*(mpc.iS_0.*mpc.Ai_0);
for k = 1:mpc.N-1
    mpc.H_k(:,:,k) = mpc.t*mpc.H_f0_k(:,:,k) + ...
                        mpc.Ai_k(:,:,k)'*(mpc.iS_k(:,k).*mpc.Ai_k(:,:,k));

    mpc.Q_k(:,:,k) = mpc.H_k(mpc.se_col,mpc.se_col,k);
    mpc.R_k(:,:,k) = mpc.H_k(mpc.u_col,mpc.u_col,k);
    mpc.Y_k(:,:,k) = mpc.H_k(mpc.u_col,mpc.se_col,k);
end

mpc.Q_ter(:,:) = mpc.t*mpc.H_f0_ter;
if mpc.ng_k(3)
    mpc.Q_ter(:,:) = mpc.Q_ter(:,:) + mpc.Ai_ter'*(mpc.iS_ter.*mpc.Ai_ter);
end

end