function [beq_0,beq_k] = update_mpc_beq(beq_0,beq_k,A,Bd,d,s_prev,dyn_use_d,N)
%% dynamics
% k = 0
beq_0(:) = -A(:,:,1)*s_prev;

if dyn_use_d
    beq_0(:) = beq_0(:) - Bd(:,:,1)*d(:,1);

    for k = 1:N-1
        beq_k(:,k) = -Bd(:,:,k+1)*d(:,k+1);
    end
end

end
