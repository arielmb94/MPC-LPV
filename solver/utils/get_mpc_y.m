% Computes output vectors of the form: y = Cs + Du + Dd
function mpc = get_mpc_y(mpc,u_k,d_k)

mpc.y(:,:)=0;
for k = 1:mpc.N-1
    if mpc.y_use_s
        mpc.y(:,k) = mpc.y(:,k) + mpc.C*mpc.s(:,k);
    end
    if mpc.y_use_u
        mpc.y(:,k) = mpc.y(:,k) + mpc.D*u_k(:,k);
    end
    if mpc.y_use_d
        mpc.y(:,k) = mpc.y(:,k) + mpc.Dd*d_k(:,k);
    end
end
mpc.err(:,:) = mpc.r-mpc.y;
end