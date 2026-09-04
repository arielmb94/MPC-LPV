function [recompute_cost_hess,R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Q_f0_ter] = ...
    update_mpc_f0_hess(recompute_cost_hess,R_f0_0,Q_f0_k,R_f0_k,Y_f0_k,Q_f0_ter,...
                        tracking_cost,quad_control_cost,controlrate_cost,...
                        quad_custom_cost,ter_ingredients,...
                        y_use_k0,y_use_s,y_use_u,y_use_ter,...
                        z_use_k0,z_use_s,z_use_su,z_use_u,z_use_ter,...
                        R_E_0,Q_E,R_E,Y_E,Q_E_ter,Ru,Rdu,...
                        R_Z_0,Q_Z,R_Z,Y_Z,Q_Z_ter,P2,s_col,su_col,N)

recompute_cost_hess = 0;

R_f0_0(:,:) = 0;
Q_f0_k(:,:,:) = 0;
R_f0_k(:,:,:) = 0;
Y_f0_k(:,:,:) = 0;
Q_f0_ter(:,:) = 0;

if tracking_cost
    if y_use_k0
        R_f0_0(:,:) = R_f0_0 + R_E_0;
    end

    if y_use_s && y_use_u
        Q_f0_k(s_col,s_col,:) = Q_f0_k(s_col,s_col,:) + Q_E;
        R_f0_k(:,:,:) = R_f0_k + R_E;
        Y_f0_k(:,s_col,:) = Y_f0_k(:,s_col,:) + Y_E;
    elseif y_use_s
        Q_f0_k(s_col,s_col,:) = Q_f0_k(s_col,s_col,:) + Q_E;
    elseif y_use_u
        R_f0_k(:,:,:) = R_f0_k + R_E;
    end

    if y_use_ter
        Q_f0_ter(s_col,s_col) = Q_f0_ter(s_col,s_col) + Q_E_ter;
    end
end

if quad_control_cost
    R_f0_0(:,:) = R_f0_0 + Ru(:,:,1);

    for k = 1:N-1
        R_f0_k(:,:,k) = R_f0_k(:,:,k) + Ru(:,:,k+1);
    end
end

if controlrate_cost
    R_f0_0(:,:) = R_f0_0 + Rdu(:,:,1);

    for k = 1:N-1
        Q_f0_k(su_col,su_col,k) = Q_f0_k(su_col,su_col,k) + Rdu(:,:,k+1);
        R_f0_k(:,:,k) = R_f0_k(:,:,k) + Rdu(:,:,k+1);
        Y_f0_k(:,su_col,k) = Y_f0_k(:,su_col,k) - Rdu(:,:,k+1);
    end
end

if quad_custom_cost
    if z_use_k0
        R_f0_0(:,:) = R_f0_0 + R_Z_0;
    end

    if z_use_s && z_use_su && z_use_u
        Q_f0_k(:,:,:) = Q_f0_k + Q_Z;
        R_f0_k(:,:,:) = R_f0_k + R_Z;
        Y_f0_k(:,:,:) = Y_f0_k + Y_Z;
    elseif z_use_s && z_use_su
        Q_f0_k(:,:,:) = Q_f0_k + Q_Z;
    elseif z_use_s && z_use_u
        Q_f0_k(s_col,s_col,:) = Q_f0_k(s_col,s_col,:) + Q_Z;
        R_f0_k(:,:,:) = R_f0_k + R_Z;
        Y_f0_k(:,s_col,:) = Y_f0_k(:,s_col,:) + Y_Z;
    elseif z_use_su && z_use_u
        Q_f0_k(su_col,su_col,:) = Q_f0_k(su_col,su_col,:) + Q_Z;
        R_f0_k(:,:,:) = R_f0_k + R_Z;
        Y_f0_k(:,su_col,:) = Y_f0_k(:,su_col,:) + Y_Z;
    elseif z_use_s
        Q_f0_k(s_col,s_col,:) = Q_f0_k(s_col,s_col,:) + Q_Z;
    elseif z_use_su
        Q_f0_k(su_col,su_col,:) = Q_f0_k(su_col,su_col,:) + Q_Z;
    elseif z_use_u
        R_f0_k(:,:,:) = R_f0_k + R_Z;
    end

    if z_use_ter
        Q_f0_ter(s_col,s_col) = Q_f0_ter(s_col,s_col) + Q_Z_ter;
    end
end

if ter_ingredients
    Q_f0_ter(s_col,s_col) = Q_f0_ter(s_col,s_col) + P2;
end

end
