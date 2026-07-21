function mpc = recover_slacks(mpc,delta_u,delta_se)
%% \mu_i = (-S)^{-1}(-\hat{r}_i-A_i\Delta x) = S^{-1}(\hat{r}_i+A_i\Delta x)

%k = 0
mpc.mu_i_0(:) = mpc.iS_0.*(mpc.ri_hat_0+mpc.Ai_0*delta_u(:,1));

for k = 1:mpc.N-1
    mpc.mu_i_k(:,k) = mpc.iS_k(:,k).*(mpc.ri_hat_k(:,k) + ... 
                                  mpc.Ai_k(:,mpc.se_col,k)*delta_se(:,k) + ...
                                  mpc.Ai_k(:,mpc.u_col,k)*delta_u(:,k+1));
end
%k = N
mpc.mu_i_ter(:) = mpc.iS_ter.*(mpc.ri_hat_ter + mpc.Ai_ter*delta_se(:,mpc.N));

%% \Delta g  = (g^2)(-r_g-\mu_i) 
%  \Delta v = (v^2)(-r_v+ \mu_i) 

mpc.delta_g_0(:) = -mpc.g2_0.*(mpc.rg_0+mpc.mu_i_0);
if mpc.nv_k(1)
    mpc.delta_v_0(:) = mpc.v2_0.*(mpc.mu_i_0(mpc.v_rows_0)-mpc.rv_0);
end

for k = 1:mpc.N-1
    mpc.delta_g_k(:,k) = -mpc.g2_k(:,k).*(mpc.rg_k(:,k)+mpc.mu_i_k(:,k));
    if mpc.nv_k(2)
        mpc.delta_v_k(:,k) = mpc.v2_k(:,k).*(mpc.mu_i_k(mpc.v_rows_k,k)-mpc.rv_k(:,k));
    end
end

if mpc.ng_k(3)
    mpc.delta_g_ter = -mpc.g2_ter.*(mpc.rg_ter+mpc.mu_i_ter);
    mpc.delta_v_ter = mpc.v2_ter.*(mpc.mu_i_ter-mpc.rv_ter);
end

end