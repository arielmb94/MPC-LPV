function mpc = get_mpc_x(mpc,x)

    mpc.s(:,:)= x(mpc.s_index_k);
    mpc.s_ter(:) = x(mpc.s_index_k(:,mpc.N));
        
end