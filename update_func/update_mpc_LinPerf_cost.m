function mpc = update_mpc_LinPerf_cost(mpc,Cz,Dz,Ddz,Qz)

mpc.recompute_cost_hess = 1;

if ~isempty(Cz)
    mpc.Cz = Cz;
end

if ~isempty(Dz)
    mpc.Dz = Dz;
end

if ~isempty(Ddz)   
    mpc.Ddz = Ddz;
end

if ~isempty(Qz)   
    mpc.Qz = Qz;
end

[mpc.gradPerfQz,mpc.hessPerfTerm] = genLinOutGradHess(mpc.Qz,mpc.Cz,mpc.Dz,...
    mpc.N,mpc.Nx,mpc.Nu,mpc.Nz,mpc.nx,mpc.nu,mpc.nz);
    
end