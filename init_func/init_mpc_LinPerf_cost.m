function mpc = set_mpc_LinPerf_cost(mpc,Qz,Cz,Dz,Ddz)

mpc.Qz = Qz;

mpc.Cz = Cz;
mpc.Dz = Dz;
mpc.Ddz = Ddz;

mpc.ndz = size(Ddz,2);  %number of disturbance inputs to performance cost
mpc.nz = size(Cz,1);  %number of performances

if mpc.Dz == 0
    mpc.Nz = (mpc.N-1)*mpc.nz;
else
    mpc.Nz = mpc.N*mpc.nz;
end

mpc.Ndz = mpc.N*mpc.ndz;

if isempty(mpc.hessCost)
    mpc.hessCost = zeros(mpc.Nu+mpc.Nx);
end

[mpc.gradPerfQz,mpc.hessPerfTerm] = genLinOutGradHess(Qz,Cz,Dz,mpc.N,...
        mpc.Nx,mpc.Nu,mpc.Nz,mpc.nx,mpc.nu,mpc.nz);
    
mpc.hessCost = mpc.hessCost + mpc.hessPerfTerm;

end