function mpc = get_mpc_diff_u(mpc,x,u_prev)

    mpc.su(:,1) = u_prev;
    mpc.su(:,2:mpc.N) = x(mpc.su_index_k(:,2:mpc.N));

    mpc.du(:,:) = mpc.u-mpc.su;
        
end