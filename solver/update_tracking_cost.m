function mpc = update_tracking_cost(mpc)

mpc.update_tracking = 0;
mpc.recompute_cost_hess = 1;

if mpc.y_use_k0
    mpc.gradErr_Qe_0(:,:) = mpc.grad_err_0*mpc.Qe_0;
    mpc.H_ErrCost_0(:,:) = mpc.gradErr_Qe_0*mpc.grad_err_0';
end

for k = 1:mpc.N-1
    mpc.gradErr_Qe_k(:,:,k) = mpc.grad_err*mpc.Qe;
    mpc.H_ErrCost_k(:,:,k) = mpc.gradErr_Qe_k(:,:,k)*mpc.grad_err';
end

if mpc.y_use_ter
    mpc.gradErr_Qe_ter(:,:) = mpc.grad_err_ter*mpc.Qe_ter;
    mpc.H_ErrCost_ter(:,:) = mpc.gradErr_Qe_ter*mpc.grad_err_ter';
end

end