function mpc = riccati_KKT(mpc,A,B,B_0,Q,Q_ter,R_0,R,Y,...
                            ru_0,ru,rs,rs_ter,rp_0,rp)

linsolve_opts = struct('SYM',true);

mpc.Q_hat_ric(:,:,mpc.N) = Q_ter;
mpc.rs_hat_ric(:,mpc.N) = rs_ter;

for k = mpc.N-1:-1:1

    % Cache Q_hat*A and Q_hat*B once for this stage.
    mpc.QA_ric(:,:) = mpc.Q_hat_ric(:,:,k+1)*A(:,:,k);
    mpc.QB_ric(:,:) = mpc.Q_hat_ric(:,:,k+1)*B(:,:,k);

    % \hat{R}_k = R_k + B^T\hat{Q}_{k+1}B
    mpc.R_hat_ric(:,:,k) = R(:,:,k) + B(:,:,k)'*mpc.QB_ric;
    for j = 1:mpc.nu
        mpc.R_hat_ric(j,j,k) = mpc.R_hat_ric(j,j,k) + mpc.eps_thknv;
    end
    % \hat{Y}_k = Y_k + B^T\hat{Q}_{k+1}A
    mpc.Y_hat_ric(:,:,k) = Y(:,:,k) + B(:,:,k)'*mpc.QA_ric;

    % \hat{r}_{pk} = \hat{r}_{sk+1} + \hat{Q}_{k+1}r_{pk}
    mpc.rp_hat_ric(:,k) = mpc.rs_hat_ric(:,k+1) + ...
                          mpc.Q_hat_ric(:,:,k+1)*rp(:,k);
    % \hat{r}_{uk} = r_{uk} + B^T\hat{r}_{pk}
    mpc.ru_hat_ric(:,k) = ru(:,k) + B(:,:,k)'*mpc.rp_hat_ric(:,k);

    % Solve [K_k, d_k] in one factorization of the regularized block.
    mpc.solve_rhs_ric(:,1:mpc.nse) = mpc.Y_hat_ric(:,:,k);
    mpc.solve_rhs_ric(:,mpc.nse+1) = mpc.ru_hat_ric(:,k);
    mpc.solve_result_ric(:,:) = linsolve(mpc.R_hat_ric(:,:,k),...
                                         mpc.solve_rhs_ric,linsolve_opts);
    mpc.K_ric(:,:,k) = mpc.solve_result_ric(:,1:mpc.nse);
    mpc.d_ric(:,k) = mpc.solve_result_ric(:,mpc.nse+1);

    % \hat{Q}_k = Q_k + A^T\hat{Q}_{k+1}A - \hat{Y}_k^TK_k
    mpc.Q_hat_ric(:,:,k) = Q(:,:,k) + A(:,:,k)'*mpc.QA_ric - ...
                           mpc.Y_hat_ric(:,:,k)'*mpc.K_ric(:,:,k);
    % \hat{r}_{sk} = r_{sk} + A^T\hat{r}_{pk} - \hat{Y}_k^Td_k
    mpc.rs_hat_ric(:,k) = rs(:,k) + A(:,:,k)'*mpc.rp_hat_ric(:,k) - ...
                          mpc.Y_hat_ric(:,:,k)'*mpc.d_ric(:,k);
end

mpc.QB_ric(:,:) = mpc.Q_hat_ric(:,:,1)*B_0;
mpc.R_hat_ric_0(:,:) = R_0+B_0'*mpc.QB_ric;
for j = 1:mpc.nu
    mpc.R_hat_ric_0(j,j) = mpc.R_hat_ric_0(j,j) + mpc.eps_thknv;
end

mpc.rp_hat_ric_0(:) = mpc.rs_hat_ric(:,1) + ...
                       mpc.Q_hat_ric(:,:,1)*rp_0;
mpc.ru_hat_ric_0(:) = ru_0 + B_0'*mpc.rp_hat_ric_0;
mpc.d_ric_0(:) = linsolve(mpc.R_hat_ric_0,mpc.ru_hat_ric_0,...
                          linsolve_opts);

% u_0 = -d_0
mpc.delta_u(:,1) = -mpc.d_ric_0;
% s_1 = Bu_0 + r_{p0}
mpc.delta_se(:,1) = B_0*mpc.delta_u(:,1) + rp_0;
% \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
%mu(:,1) = mpc.rs_hat_ric(:,1) + mpc.Q_hat_ric(:,:,1)*mpc.delta_se(:,1);

for k = 1:mpc.N-1
    ku = k+1; % u mu vectors shifted by 1 stage
    % u_k = -d_k - K_ks_k
    mpc.delta_u(:,ku) = -mpc.d_ric(:,k) - ...
                         mpc.K_ric(:,:,k)*mpc.delta_se(:,k);
    % s_{k+1} =  As_k + Bu_k + r_{pk}
    mpc.delta_se(:,k+1) = A(:,:,k)*mpc.delta_se(:,k) + B(:,:,k)*mpc.delta_u(:,ku) + rp(:,k);
    % \mu_k = \hat{r}_{sk+1} + \hat{Q}_{k+1}s_{k+1}
    %mu(:,ku) = mpc.rs_hat_ric(:,k+1) + mpc.Q_hat_ric(:,:,k+1)*mpc.delta_se(:,k+1);
end

end

