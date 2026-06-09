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
   mpc.slacks(:) = x(mpc.slack_index);
   mpc.g(:) = x(mpc.g_index);
end
if any(mpc.nv_k)
   mpc.v(:) = x(mpc.v_index);
end

end