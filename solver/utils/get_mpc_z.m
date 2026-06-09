% Computes vectors: h = C*s + D*u Dsu*su + Dd*d
function mpc = get_mpc_z(mpc)

mpc.z(:,:)=0;

%k=0 (only if z depends on u)
if mpc.z_use_u
    mpc.z(:,1) = mpc.Dz*mpc.u(:,1);
    if mpc.z_use_s
        mpc.z(:,1) = mpc.z(:,1) + mpc.Cz*mpc.s(:,1);
    end
    if mpc.z_use_su
        mpc.z(:,1) = mpc.z(:,1) + mpc.Dsuz*mpc.su(:,1);
    end
    if mpc.z_use_d
        mpc.z(:,1) = mpc.z(:,1) + mpc.Ddz*mpc.dz(:,1);
    end
end

for k = 2:mpc.N
    if mpc.z_use_s
        mpc.z(:,k) = mpc.z(:,k) + mpc.Cz*mpc.s(:,k);
    end
    if mpc.z_use_u
        mpc.z(:,k) = mpc.z(:,k) + mpc.Dz*mpc.u(:,k);
    end
    if mpc.z_use_su
        mpc.z(:,k) = mpc.z(:,k) + mpc.Dsuz*mpc.su(:,k);
    end
    if mpc.z_use_d
        mpc.z(:,k) = mpc.z(:,k) + mpc.Ddz*mpc.dz(:,k);
    end
end
end