function mpc = init_mpc_DiffControl_cost(mpc,Rdu)

mpc.Rdu = Rdu;

if isempty(mpc.hessCost)
    mpc.hessCost = zeros(mpc.Nu+mpc.Nx);
end

[mpc.gradDiffCtlrRdu,mpc.hessDiffCtrlTerm] = genDiffControlGradHess(Rdu,mpc.N,...
    mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

mpc.hessCost = mpc.hessCost + mpc.hessDiffCtrlTerm;

end