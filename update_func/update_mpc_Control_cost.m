function mpc = update_mpc_Control_cost(mpc,Ru)

mpc.Ru = Ru;

[mpc.gradCtlrRu,mpc.hessCtrlTerm] = genControlGradHess(Ru,mpc.N,...
    mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

mpc.recompute_cost_hess = 1;

end