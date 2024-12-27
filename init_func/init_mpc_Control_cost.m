function mpc = init_mpc_Control_cost(mpc,Ru)

mpc.Ru = Ru;

if isempty(mpc.hessCost)
    mpc.hessCost = zeros(mpc.Nu+mpc.Nx);
end

[mpc.gradCtlrRu,mpc.hessCtrlTerm] = genControlGradHess(Ru,mpc.N_ctr_hor,...
    mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

mpc.hessCost = mpc.hessCost + mpc.hessCtrlTerm;

end