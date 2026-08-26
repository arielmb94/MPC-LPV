% gradients exclusively on the cost function on s, su and u (no slacks)
function [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = grad_f0_MPC(mpc,...
                           grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,...
                           grad_se_f0_ter,s_col,su_col,N,t)


grad_u_f0_0(:) = 0;
grad_se_f0_k(:,:) = 0;
grad_u_f0_k(:,:) = 0;
grad_se_f0_ter(:) = 0;

if mpc.quad_control_cost || mpc.lin_control_cost
    [grad_u_f0_0,grad_u_f0_k] = control_f0(mpc.quad_control_cost,mpc.lin_control_cost,...
                                        grad_u_f0_0,grad_u_f0_k,mpc.Ru,mpc.ru,mpc.u,N);
end

if mpc.controlrate_cost
    [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k] = controlrate_f0(grad_u_f0_0,...
                                    grad_se_f0_k,grad_u_f0_k,...
                                    mpc.Rdu,mpc.du,su_col,N);
end

if mpc.tracking_cost
    [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = tracking_f0(...
                    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
                    mpc.y_use_k0,mpc.y_use_ter,mpc.y_use_s,mpc.y_use_u,...
                    mpc.grad_u_E_0,mpc.grad_s_E,mpc.grad_u_E,mpc.grad_s_E_ter,...
                    mpc.err_0,mpc.err,mpc.err_ter,N,s_col);
end

if mpc.quad_custom_cost
    [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = custom_quad_f0(...
                    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
                    mpc.z_use_k0,mpc.z_use_ter,mpc.z_use_s,mpc.z_use_su,mpc.z_use_u,...
                    mpc.grad_u_Z_0,mpc.grad_s_Z,mpc.grad_su_Z,mpc.grad_u_Z,mpc.grad_s_Z_ter,...
                    mpc.z_0,mpc.z,mpc.z_ter,N,s_col,su_col);

end

if mpc.lin_custom_cost
    [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = custom_lin_f0(...
                    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
                    mpc.z_use_k0,mpc.z_use_ter,mpc.z_use_s,mpc.z_use_su,mpc.z_use_u,...
                    mpc.grad_u_Zlin_0,mpc.grad_s_Zlin,mpc.grad_su_Zlin,...
                    mpc.grad_u_Zlin,mpc.grad_s_Zlin_ter,s_col,su_col);

end

if mpc.ter_ingredients
    grad_se_f0_ter(s_col) = grad_se_f0_ter(s_col) - mpc.P2*(mpc.xN_ref-mpc.s_ter);
end

grad_u_f0_0 = t*grad_u_f0_0;
grad_se_f0_k = t*grad_se_f0_k;
grad_u_f0_k = t*grad_u_f0_k;
grad_se_f0_ter = t*grad_se_f0_ter;
end




function [grad_u_f0_0,grad_u_f0_k] = control_f0(quad_control_cost,lin_control_cost,...
                                    grad_u_f0_0,grad_u_f0_k,Ru,ru,u,N)
    if quad_control_cost
        grad_u_f0_0 = grad_u_f0_0 + Ru(:,:,1)*u(:,1);
    
        for k = 1:N-1
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + Ru(:,:,k+1)*u(:,k+1);
        end
    end
    if lin_control_cost
        grad_u_f0_0 = grad_u_f0_0 + ru(:,1);
    
        grad_u_f0_k = grad_u_f0_k + ru(:,2:N);
    end
end

function [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k] = controlrate_f0(grad_u_f0_0,...
                                    grad_se_f0_k,grad_u_f0_k,Rdu,du,su_idx,N)
    % k = 0
    grad_u_f0_0 = grad_u_f0_0 + Rdu(:,:,1)*du(:,1);

    for k = 1:N-1
        grad_du_k = Rdu(:,:,k+1)*du(:,k+1);

        grad_se_f0_k(su_idx,k) = grad_se_f0_k(su_idx,k) - grad_du_k;
        grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_du_k;
    end
end

function [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = tracking_f0(...
                    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
                    y_use_k0,y_use_ter,y_use_s,y_use_u,...
                    grad_u_E_0,grad_s_E_k,grad_u_E_k,grad_s_E_ter,...
                    err_0,err,err_ter,N,s_idx)

    if y_use_k0
        grad_u_f0_0 = grad_u_f0_0 - grad_u_E_0*err_0;
    end

    for k = 1:N-1
        if y_use_s && y_use_u
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) -...
                                        grad_s_E_k(:,:,k)*err(:,k);
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) -...
                                        grad_u_E_k(:,:,k)*err(:,k);
        elseif y_use_s
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) -...
                                        grad_s_E_k(:,:,k)*err(:,k);
        elseif y_use_u
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) -...
                                        grad_u_E_k(:,:,k)*err(:,k);
        end
    end

    if y_use_ter
        grad_se_f0_ter(s_idx) = grad_se_f0_ter(s_idx) - grad_s_E_ter*err_ter;
    end
end

function [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = custom_quad_f0(...
                    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
                    z_use_k0,z_use_ter,z_use_s,z_use_su,z_use_u,...
                    grad_u_Z_0,grad_s_Z,grad_su_Z,grad_u_Z,grad_s_Z_ter,...
                    z_0,z,z_ter,N,s_idx,su_idx)

    if z_use_k0
        grad_u_f0_0 = grad_u_f0_0 + grad_u_Z_0*z_0;
    end

    for k = 1:N-1
        if z_use_s && z_use_su && z_use_u

            z_k = z(:,k);
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) + grad_s_Z(:,:,k)*z_k;
            grad_se_f0_k(su_idx,k) = grad_se_f0_k(su_idx,k) + grad_su_Z(:,:,k)*z_k;
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z_k;
        elseif z_use_s && z_use_su

            z_k = z(:,k);
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) + grad_s_Z(:,:,k)*z_k;
            grad_se_f0_k(su_idx,k) = grad_se_f0_k(su_idx,k) + grad_su_Z(:,:,k)*z_k;
        elseif z_use_s && z_use_u

            z_k = z(:,k);
            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) + grad_s_Z(:,:,k)*z_k;
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z_k;
        elseif z_use_su && z_use_u

            z_k = z(:,k);
            grad_se_f0_k(su_idx,k) = grad_se_f0_k(su_idx,k) + grad_su_Z(:,:,k)*z_k;
            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z_k;
        elseif z_use_s

            grad_se_f0_k(s_idx,k) = grad_se_f0_k(s_idx,k) + grad_s_Z(:,:,k)*z(:,k);
        elseif z_use_su

            grad_se_f0_k(su_idx,k) = grad_se_f0_k(su_idx,k) + grad_su_Z(:,:,k)*z(:,k);
        elseif z_use_u

            grad_u_f0_k(:,k) = grad_u_f0_k(:,k) + grad_u_Z(:,:,k)*z(:,k);
        end
    end

    if z_use_ter
        grad_se_f0_ter(s_idx) = grad_se_f0_ter(s_idx) + grad_s_Z_ter*z_ter;
    end
end

function [grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter] = custom_lin_f0(...
                    grad_u_f0_0,grad_se_f0_k,grad_u_f0_k,grad_se_f0_ter,...
                    z_use_k0,z_use_ter,z_use_s,z_use_su,z_use_u,...
                    grad_u_Z_0,grad_s_Z,grad_su_Z,grad_u_Z,grad_s_Z_ter,...
                    s_idx,su_idx)

if z_use_k0
    grad_u_f0_0 = grad_u_f0_0 + grad_u_Z_0;
end

% k = 1:N-1
if z_use_s && z_use_su && z_use_u

    grad_se_f0_k(s_idx,:) = grad_se_f0_k(s_idx,:) + grad_s_Z;
    grad_se_f0_k(su_idx,:) = grad_se_f0_k(su_idx,:) + grad_su_Z;
    grad_u_f0_k = grad_u_f0_k + grad_u_Z;
elseif z_use_s && z_use_su

    grad_se_f0_k(s_idx,:) = grad_se_f0_k(s_idx,:) + grad_s_Z;
    grad_se_f0_k(su_idx,:) = grad_se_f0_k(su_idx,:) + grad_su_Z;
elseif z_use_s && z_use_u

    grad_se_f0_k(s_idx,:) = grad_se_f0_k(s_idx,:) + grad_s_Z;
    grad_u_f0_k = grad_u_f0_k + grad_u_Z;
elseif z_use_su && z_use_u

    grad_se_f0_k(su_idx,:) = grad_se_f0_k(su_idx,:) + grad_su_Z;
    grad_u_f0_k = grad_u_f0_k + grad_u_Z;
elseif z_use_s

    grad_se_f0_k(s_idx,:) = grad_se_f0_k(s_idx,:) + grad_s_Z;
elseif z_use_su

    grad_se_f0_k(su_idx,:) = grad_se_f0_k(su_idx,:) + grad_su_Z;
elseif z_use_u

    grad_u_f0_k = grad_u_f0_k + grad_u_Z;
end

if z_use_ter
    grad_se_f0_ter(s_idx) = grad_se_f0_ter(s_idx) + grad_s_Z_ter;
end

end


