% Computes output vectors of the form: y = Cs + Du + Dd
function mpc = get_mpc_y(mpc)

mpc.y(:,:)=0;
for k = 2:mpc.N
    if mpc.y_use_s
        mpc.y(:,k-1) = mpc.y(:,k-1) + mpc.C*mpc.s(:,k);
    end
    if mpc.y_use_u
        mpc.y(:,k-1) = mpc.y(:,k-1) + mpc.D*mpc.u(:,k);
    end
    if mpc.y_use_d
        mpc.y(:,k-1) = mpc.y(:,k-1) + mpc.Dd*mpc.d(:,k);
    end
end
mpc.err = mpc.y-mpc.r;
end