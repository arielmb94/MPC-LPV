function [delta_x,grad_J] = stage2vec(d_u,d_se,d_g_0,d_g_k,d_g_ter,...
                                      d_v_0,d_v_k,d_v_ter,...
                                      g_0,g_k,g_ter,...
                                      v_0,v_k,v_ter,...
                                      t,grad_qv_0,grad_qv_k,grad_qv_ter,...
                                      grad_u_f0_0,grad_u_f0_k,...
                                      grad_se_f0_k,grad_se_f0_ter,...
                                      u_index_k,se_index_k,...
                                      g_index_0,g_index_k,g_index_ter,...
                                      v_index_0,v_index_k,v_index_ter,...
                                      ng_k,nv_k,n,N)

delta_x = zeros(n,1);

delta_x(u_index_k) = d_u;
delta_x(se_index_k) = d_se;

delta_x(g_index_0) = d_g_0;
delta_x(g_index_k) = d_g_k;
delta_x(g_index_ter) = d_g_ter;

delta_x(v_index_0) = d_v_0;
delta_x(v_index_k) = d_v_k;
delta_x(v_index_ter) = d_v_ter;


grad_J = zeros(1,n);

grad_J(u_index_k(:,1)) = grad_u_f0_0;
grad_J(u_index_k(:,2:N)) = grad_u_f0_k;

grad_J(se_index_k(:,1:N-1)) = grad_se_f0_k;
grad_J(se_index_k(:,N)) = grad_se_f0_ter;

if ng_k(1), grad_J(g_index_0) = -1./g_0; end
if ng_k(2), grad_J(g_index_k) = -1./g_k; end
if ng_k(3), grad_J(g_index_ter) = -1./g_ter; end

if nv_k(1), grad_J(v_index_0) = t*grad_qv_0-1./v_0; end
if nv_k(2), grad_J(v_index_k) = t*grad_qv_k-1./v_k; end
if nv_k(3), grad_J(v_index_ter) = t*grad_qv_ter-1./v_ter; end

end