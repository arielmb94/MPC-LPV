function mpc = get_mpc_u(mpc,x)

    mpc.u(:,:) = x(mpc.u_index_k);
        
end