function [delta_var,delta_g,delta_v] = reduced_KKT(mpc,x0,grad_f0,opts)

hess_f0 = mpc.t*mpc.hessCost(mpc.variables_index,mpc.variables_index);
grad_f0_x = mpc.t*grad_f0(mpc.variables_index); 

Adyn = [mpc.Aeq(mpc.index_eq_k,mpc.variables_index);...
       mpc.Aeq(mpc.index_eq_N,mpc.variables_index)];
Ainq = [mpc.Aeq(mpc.index_ineq_0,mpc.variables_index);...
       mpc.Aeq(mpc.index_ineq_k,mpc.variables_index);...
       mpc.Aeq(mpc.index_ineq_ter,mpc.variables_index)];

rp = mpc.Aeq*x0-mpc.beq;

rp_dyn = [];
rp_inq = [];

for k = 1:mpc.N-1
    rp_dyn = [rp_dyn;rp(mpc.index_eq_k(:,k))];
end
rp_dyn = [rp_dyn;rp(mpc.index_eq_N)];

rp_inq = [rp_inq;rp(mpc.index_ineq_0)];
for k = 1:mpc.N-1
    rp_inq = [rp_inq;rp(mpc.index_ineq_k(:,k))];
end
rp_inq = [rp_inq;rp(mpc.index_ineq_ter)];

m = sum(mpc.ng_k);
S = zeros(m,1);
rg = [];
rv = [];
ri_hat = zeros(m,1);

ineq_index = cumsum(mpc.ng_k);

% form S, ri hat
% ri_hat = ri -g^2rg+v^2rv

%k =0
index_ineq_0 = 1:ineq_index(1);

g0 = x0(mpc.S_gi_0)-mpc.slack_epsilon;
rg0 = -1./g0;
rg = [rg;rg0];
S(index_ineq_0) = g0.^2;

ri_hat(index_ineq_0) = rp_inq(index_ineq_0) - (g0.^2).*rg0;

if any(mpc.S_vi_0)
    ng_k = mpc.ng_k(1);
    vk = zeros(ng_k,1);
    rvk = zeros(ng_k,1);

    vk(mpc.S_vi_0_tmplt) = x0(mpc.S_vi_0)-mpc.slack_epsilon;
    rvk(mpc.S_vi_0_tmplt) = mpc.t*grad_f0(mpc.S_vi_0) - 1./vk(mpc.S_vi_0_tmplt);
    rv = [rv;rvk(mpc.S_vi_0_tmplt)];

    S(index_ineq_0) = S(index_ineq_0) + vk.^2;

    ri_hat(index_ineq_0) = ri_hat(index_ineq_0) + (vk.^2).*rvk;
end

for k = 1:mpc.N-1
    index_ineq_k = ineq_index(k)+1:ineq_index(k+1);
    gk = x0(mpc.S_gi_k(:,k))-mpc.slack_epsilon;
    S(index_ineq_k) = gk.^2;

    rgk = -1./gk;
    rg = [rg;rgk];
    ri_hat(index_ineq_k) = rp_inq(index_ineq_k) - (gk.^2).*rgk;

    if any(mpc.S_vi_k(:,k))
        ng_k = mpc.ng_k(k+1);
        vk = zeros(ng_k,1);
        rvk = zeros(ng_k,1);

        % crear: index_ineq_k(mpc.S_vi_k_tmplt)

        vk(mpc.S_vi_k_tmplt) = x0(mpc.S_vi_k(:,k))-mpc.slack_epsilon;
        rvk(mpc.S_vi_k_tmplt) = mpc.t*grad_f0(mpc.S_vi_k(:,k)) - 1./vk(mpc.S_vi_k_tmplt);
        rv = [rv;rvk(mpc.S_vi_k_tmplt)];

        S(index_ineq_k) = S(index_ineq_k)+ vk.^2;

        ri_hat(index_ineq_k) = ri_hat(index_ineq_k) + vk.^2.*rvk;
    end
end

if any(mpc.S_vi_ter)
    index_ineq_ter = ineq_index(mpc.N)+1:ineq_index(mpc.N+1);
    gk = x0(mpc.S_gi_ter)-mpc.slack_epsilon;
    rgk = -1./gk;
    rg = [rg;rgk];
    S(index_ineq_ter) = gk.^2;

    vk = x0(mpc.S_vi_ter)-mpc.slack_epsilon;
    rvk = mpc.t*grad_f0(mpc.S_vi_ter) - 1./vk;
    rv = [rv;rvk];
    S(index_ineq_ter) = S(index_ineq_ter)+ vk.^2;

    ri_hat(index_ineq_ter) = rp_inq(index_ineq_ter) - gk.^2.*rgk + vk.^2.*rvk;
end

iS = 1./S;

rx_hat = grad_f0_x + Ainq'*diag(iS)*ri_hat;


neq = size(Adyn,1);
nvar = size(Adyn,2);

H = hess_f0 + Ainq'*diag(iS)*Ainq+eye(nvar)*mpc.eps_thknv;


KKT = [H Adyn';
       Adyn zeros(neq,neq)];

delta_x = - linsolve(KKT,[rx_hat;rp_dyn],opts);

% recover solution

delta_var = delta_x(1:nvar);

mu_i = (ri_hat+Ainq*delta_var)./S;
delta_g = ((mpc.g-mpc.slack_epsilon).^2).*(1./(mpc.g-mpc.slack_epsilon)-mu_i);

mu_vi = [];
if any(mpc.S_vi_0)
    index_ineq_0 = 1:ineq_index(1);
    mu_vi_k = mu_i(index_ineq_0(mpc.S_vi_0_tmplt));
    mu_vi = [mu_vi;mu_vi_k];
end
for k = 1:mpc.N-1
    index_ineq_k = ineq_index(k)+1:ineq_index(k+1);
    mu_vi_k = mu_i(index_ineq_k(mpc.S_vi_k_tmplt));
    mu_vi = [mu_vi;mu_vi_k];
end
if any(mpc.S_vi_ter)
    index_ineq_ter = ineq_index(mpc.N)+1:ineq_index(mpc.N+1);
    mu_vi_k = mu_i(index_ineq_ter);
    mu_vi = [mu_vi;mu_vi_k];
end

delta_v = (mpc.v-mpc.slack_epsilon).^2.*(-rv+mu_vi);

end