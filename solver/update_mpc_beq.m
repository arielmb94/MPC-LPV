function mpc = update_mpc_beq(mpc,s_prev)
%% dynamics
% k = 0
mpc.beq_0(:) = -mpc.A(:,:,1)*s_prev;

if mpc.dyn_use_d
    mpc.beq_0(:) = mpc.beq_0(:) - mpc.Bd(:,:,1)*mpc.d(:,1);

    for k = 1:mpc.N-1
        mpc.beq_k(:,k) = - mpc.Bd(:,:,k+1)*mpc.d(:,k+1);
    end
end
