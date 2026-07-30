% Computes output vectors of the form: y = Cs + Du + Dd
function mpc = get_mpc_y(mpc,s_prev)

if mpc.y_use_k0
    mpc.y_0(:) = mpc.D_0*mpc.u(:,1);
    if mpc.y_use_s
        mpc.y_0(:) = mpc.y_0(:) + mpc.C_0*s_prev;
    end
    if mpc.y_use_d
        mpc.y_0(:) = mpc.y_0(:) + mpc.Dd_0*mpc.d(:,1);
    end

    mpc.err_0(:) = mpc.r_0 - mpc.y_0;
end

mpc.y(:,:)=0;
for k = 1:mpc.N-1
    if mpc.y_use_s
        mpc.y(:,k) = mpc.y(:,k) + mpc.C*mpc.s(:,k);
    end
    if mpc.y_use_u
        mpc.y(:,k) = mpc.y(:,k) + mpc.D*mpc.u(:,k+1);
    end
    if mpc.y_use_d
        mpc.y(:,k) = mpc.y(:,k) + mpc.Dd*mpc.d(:,k+1);
    end
end
mpc.err(:,:) = mpc.r-mpc.y;

if mpc.y_use_ter
    mpc.y_ter(:) = mpc.C_ter*mpc.s(:,mpc.N);

    mpc.err_ter(:) = mpc.r_ter - mpc.y_ter;
end
end