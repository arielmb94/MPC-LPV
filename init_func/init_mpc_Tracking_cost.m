function mpc = set_mpc_Tracking_cost(mpc,Qe)

mpc.Qe = Qe;

if isempty(mpc.hessCost)
    mpc.hessCost = zeros(mpc.Nu+mpc.Nx);
end

[mpc.gradErrQe,mpc.hessErrTerm] = genLinOutGradHess(Qe,mpc.C,mpc.D,mpc.N,...
        mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);

mpc.hessCost = mpc.hessCost + mpc.hessErrTerm;

end