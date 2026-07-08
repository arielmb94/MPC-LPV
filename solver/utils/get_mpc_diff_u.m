function mpc = get_mpc_diff_u(mpc,x,u_prev)

    mpc.su(:,1) = u_prev;
    mpc.su(:,2:mpc.N+1) = x(mpc.su_index_k);

    mpc.du(:,:) = mpc.u-mpc.su(:,1:mpc.N);
        
end