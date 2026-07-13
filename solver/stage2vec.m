function [delta_x,grad_J] = stage2vec(mpc,u,se,g_0,g_k,g_ter,v_0,v_k,v_ter)

delta_x = zeros(mpc.n,1);

delta_x(mpc.u_index_k) = u;
delta_x(mpc.se_index_k) = se;

delta_x(mpc.g_index_0) = g_0;
delta_x(mpc.g_index_k) = g_k;
delta_x(mpc.g_index_ter) = g_ter;

delta_x(mpc.v_index_0) = v_0;
delta_x(mpc.v_index_k) = v_k;
delta_x(mpc.v_index_ter) = v_ter;


grad_J = zeros(1,mpc.n);

grad_J(mpc.u_index_k(:,1)) = mpc.ru_0;
grad_J(mpc.u_index_k(:,2:mpc.N)) = mpc.ru_k;

grad_J(mpc.se_index_k(:,1:mpc.N-1)) = mpc.rse_k;
grad_J(mpc.se_index_k(:,mpc.N)) = mpc.rse_ter;

grad_J(mpc.g_index_0) = mpc.rg_0;
grad_J(mpc.g_index_k) = mpc.rg_k;
grad_J(mpc.g_index_ter) = mpc.rg_ter;

grad_J(mpc.v_index_0) = mpc.rv_0;
grad_J(mpc.v_index_k) = mpc.rv_k;
grad_J(mpc.v_index_ter) = mpc.rv_ter;


end