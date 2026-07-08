function mpc = get_mpc_x(mpc,x,s_prev)

    mpc.s(:,1) = s_prev;
    mpc.s(:,2:mpc.N+1)= x(mpc.s_index_k);
    mpc.s_ter(:) = mpc.s(:,mpc.N+1);
        
end