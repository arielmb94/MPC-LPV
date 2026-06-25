function mpc = get_mpc_variables(mpc,x,s_prev,u_prev)

%states
mpc = get_mpc_x(mpc,x,s_prev);

% control actions
mpc = get_mpc_u(mpc,x);

% differential control action
if mpc.has_du
    mpc = get_mpc_diff_u(mpc,x,u_prev);
end

% tracking outputs
if mpc.ny
    mpc = get_mpc_y(mpc);
end

% general constraints
if mpc.has_h_cnstr
    mpc = get_mpc_h(mpc);
end

if mpc.quad_custom_cost || mpc.lin_custom_cost
    % compute vector z
    mpc = get_mpc_z(mpc);
end

% get slack variables 
if any(mpc.ng_k)
   %mpc.slacks(:) = x(mpc.slack_index);
   mpc.g_0(:) = x(mpc.g_index_0);
   mpc.g_k(:) = x(mpc.g_index_k);
   mpc.g_ter(:) = x(mpc.g_index_ter);
end
if any(mpc.nv_k)
   mpc.v_0(:) = x(mpc.v_index_0);
   mpc.v_k(:) = x(mpc.v_index_k);
   mpc.v_ter(:) = x(mpc.v_index_ter);
end

end