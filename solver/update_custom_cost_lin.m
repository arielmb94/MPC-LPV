function mpc = update_custom_cost_lin(mpc)

mpc.update_customcost_lin = 0;

if mpc.z_use_k0
    mpc.gradz_qz_0(:) = mpc.grad_z_0*mpc.qz_0;
end

for k = 1:mpc.N-1
    mpc.gradz_qz_k(:,k) = mpc.grad_z(:,:,k)*mpc.qz(:,k);
end

if mpc.z_use_ter
    mpc.gradz_qz_ter(:) = mpc.grad_z_ter*mpc.qz_ter;
end

end