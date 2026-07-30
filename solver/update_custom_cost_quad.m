function mpc = update_custom_cost_quad(mpc)

mpc.update_customcost_quad = 0;
mpc.recompute_cost_hess = 1;

if mpc.z_use_k0
    mpc.gradz_Qz_0(:,:) = mpc.grad_z_0*mpc.Qz_0;
    mpc.H_CustomCost_0(:,:) = mpc.gradz_Qz_0*mpc.grad_z_0';
end

for k = 1:mpc.N-1
    mpc.gradz_Qz_k(:,:,k) = mpc.grad_z*mpc.Qz;
    mpc.H_CustomCost_k(:,:,k) = mpc.gradz_Qz_k(:,:,k)*mpc.grad_z';
end

if mpc.z_use_ter
    mpc.gradz_Qz_ter(:,:) = mpc.grad_z_ter*mpc.Qz_ter;
    mpc.H_CustomCost_ter(:,:) = mpc.gradz_Qz_ter*mpc.grad_z_ter';
end

end