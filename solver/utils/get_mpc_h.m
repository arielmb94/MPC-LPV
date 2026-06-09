% Computes vectors: h = C*s + D*u Dsu*su + Dd*d
function mpc = get_mpc_h(mpc)

mpc.h(:,:)=0;

%k=0 (only if h depends on u)
if mpc.h_cnstr.use_u
    mpc.h(:,1) = mpc.Dh*mpc.u(:,1);
    if mpc.h_cnstr.use_s
        mpc.h(:,1) = mpc.h(:,1) + mpc.Ch*mpc.s(:,1);
    end
    if mpc.h_cnstr.use_su
        mpc.h(:,1) = mpc.h(:,1) + mpc.Dsuh*mpc.su(:,1);
    end
    if mpc.h_cnstr.use_d
        mpc.h(:,1) = mpc.h(:,1) + mpc.Ddh*mpc.dh(:,1);
    end
end

for k = 2:mpc.N
    if mpc.h_cnstr.use_s
        mpc.h(:,k) = mpc.h(:,k) + mpc.Ch*mpc.s(:,k);
    end
    if mpc.h_cnstr.use_u
        mpc.h(:,k) = mpc.h(:,k) + mpc.Dh*mpc.u(:,k);
    end
    if mpc.h_cnstr.use_su
        mpc.h(:,k) = mpc.h(:,k) + mpc.Dsuh*mpc.su(:,k);
    end
    if mpc.h_cnstr.use_d
        mpc.h(:,k) = mpc.h(:,k) + mpc.Ddh*mpc.dh(:,k);
    end
end
end