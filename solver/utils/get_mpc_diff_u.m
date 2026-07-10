function mpc = get_mpc_diff_u(mpc,x,u_prev)

    mpc.su(:,:) = x(mpc.su_index_k);

    mpc.du(:,1) = mpc.u(:,1)-u_prev;
    mpc.du(:,2:mpc.N) = mpc.u(:,2:mpc.N)-mpc.u(:,1:mpc.N-1);
        
end