function [delta_g_0,delta_g_k,delta_g_ter,delta_v_0,delta_v_k,delta_v_ter] = recover_slacks(mpc,delta_u,delta_se)

se_col = [mpc.inq_s_col mpc.inq_su_col];
u_col = mpc.inq_u_col;
%% \mu_i = (-S)^{-1}(-\hat{r}_i-A_i\Delta x) = S^{-1}(\hat{r}_i+A_i\Delta x)

%k = 0
mu_i_0 = mpc.iS_0.*(mpc.ri_hat_0+mpc.Ai_0*delta_u(:,1));

for k = 1:mpc.N-1
    mu_i_k(:,k) = mpc.iS_k(:,k).*(mpc.ri_hat_k(:,k) + ... 
                                  mpc.Ai_k(:,se_col,k)*delta_se(:,k) + ...
                                  mpc.Ai_k(:,u_col,k)*delta_u(:,k+1));
end
%k = N
mu_i_ter = mpc.iS_ter.*(mpc.ri_hat_ter + mpc.Ai_ter*delta_se(:,mpc.N));

%% \Delta g  = (g^2)(-r_g-\mu_i) 
%  \Delta v = (v^2)(-r_v+ \mu_i) 

delta_g_0 = -mpc.g2_0.*(mpc.rg_0+mu_i_0);
if mpc.nv_k(1)
    delta_v_0 = mpc.v2_0.*(mu_i_0(mpc.vi_0)-mpc.rv_0);
else
    delta_v_0 = [];
end

for k = 1:mpc.N-1
    delta_g_k(:,k) = -mpc.g2_k(:,k).*(mpc.rg_k(:,k)+mu_i_k(:,k));
    if mpc.nv_k(2)
        delta_v_k(:,k) = mpc.v2_k(:,k).*(mu_i_k(mpc.vi_k,k)-mpc.rv_k(:,k));
    end
end

if mpc.has_s_cnstr
    delta_g_ter = -mpc.g2_ter.*(mpc.rg_ter+mu_i_ter);
    delta_v_ter = mpc.v2_ter.*(mu_i_ter-mpc.rv_ter);
end

end