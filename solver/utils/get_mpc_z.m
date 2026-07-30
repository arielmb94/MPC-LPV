% Computes vectors: h = C*s + D*u Dsu*su + Dd*d
function mpc = get_mpc_z(mpc,s_prev,u_prev)


%k=0 (only if z depends on u)
if mpc.z_use_k0
    mpc.z_0(:) = mpc.Dz_0*mpc.u(:,1);
    if mpc.z_use_s
        mpc.z_0(:) = mpc.z_0 + mpc.Cz_0*s_prev;
    end
    if mpc.z_use_su
        mpc.z_0(:) = mpc.z_0 + mpc.Dsuz_0*u_prev;
    end
    if mpc.z_use_d
        mpc.z_0(:) = mpc.z_0 + mpc.Ddz_0*mpc.dz(:,1);
    end
end

mpc.z(:,:)=0;
for k = 1:mpc.N-1
    if mpc.z_use_s
        mpc.z(:,k) = mpc.z(:,k) + mpc.Cz*mpc.s(:,k);
    end
    if mpc.z_use_u
        mpc.z(:,k) = mpc.z(:,k) + mpc.Dz*mpc.u(:,k+1);
    end
    if mpc.z_use_su
        mpc.z(:,k) = mpc.z(:,k) + mpc.Dsuz*mpc.su(:,k);
    end
    if mpc.z_use_d
        mpc.z(:,k) = mpc.z(:,k) + mpc.Ddz*mpc.dz(:,k+1);
    end
end

if mpc.z_use_ter
    mpc.z_ter(:) = mpc.Cz_ter*mpc.s(:,mpc.N);
end

end