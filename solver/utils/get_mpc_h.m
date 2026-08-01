% Computes vectors: h = C*s + D*u Dsu*su + Dd*d
function mpc = get_mpc_h(mpc,s_prev,u_prev)

%k=0 (only if h depends on u)
if mpc.h_cnstr.use_k0
    mpc.h_0(:) = mpc.Dh_0*mpc.u(:,1);
    if mpc.h_cnstr.use_s
        mpc.h_0(:) = mpc.h_0(:) + mpc.Ch_0*s_prev;
    end
    if mpc.h_cnstr.use_su
        mpc.h_0(:) = mpc.h_0(:) + mpc.Dsuh_0*u_prev;
    end
    if mpc.h_cnstr.use_d
        mpc.h_0(:) = mpc.h_0(:) + mpc.Ddh_0*mpc.dh(:,1);
    end
end

mpc.h(:,:)=0;
for k = 1:mpc.N-1

    if mpc.h_cnstr.use_s
        mpc.h(:,k) = mpc.h(:,k) + mpc.Ch(:,:,k)*mpc.s(:,k);
    end
    if mpc.h_cnstr.use_u
        mpc.h(:,k) = mpc.h(:,k) + mpc.Dh(:,:,k)*mpc.u(:,k+1);
    end
    if mpc.h_cnstr.use_su
        mpc.h(:,k) = mpc.h(:,k) + mpc.Dsuh(:,:,k)*mpc.su(:,k);
    end
    if mpc.h_cnstr.use_d
        mpc.h(:,k) = mpc.h(:,k) + mpc.Ddh(:,:,k)*mpc.dh(:,k+1);
    end
end

if mpc.h_cnstr.use_ter
    mpc.h_ter(:) = mpc.Ch_ter*mpc.s(:,mpc.N);
end
end