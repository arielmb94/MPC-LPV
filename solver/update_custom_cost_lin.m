function [update_customcost_lin,...
          grad_u_Zlin_0,grad_s_Zlin,grad_su_Zlin,...
          grad_u_Zlin,grad_s_Zlin_ter] = update_custom_cost_lin(...
          update_customcost_lin,...
          grad_u_Zlin_0,grad_s_Zlin,grad_su_Zlin,...
          grad_u_Zlin,grad_s_Zlin_ter,...
          z_use_k0,z_use_s,z_use_su,z_use_u,z_use_ter,...
          qz_0,qz,qz_ter,Dz_0,Cz,Dsuz,Dz,Cz_ter,N)

update_customcost_lin = 0;

if z_use_k0
    grad_u_Zlin_0(:) = Dz_0'*qz_0;
end

if z_use_s && z_use_su && z_use_u
    for k = 1:N-1
        grad_s_Zlin(:,k) = Cz(:,:,k)'*qz(:,k);
        grad_su_Zlin(:,k) = Dsuz(:,:,k)'*qz(:,k);
        grad_u_Zlin(:,k) = Dz(:,:,k)'*qz(:,k);
    end
elseif z_use_s && z_use_su
    for k = 1:N-1
        grad_s_Zlin(:,k) = Cz(:,:,k)'*qz(:,k);
        grad_su_Zlin(:,k) = Dsuz(:,:,k)'*qz(:,k);
    end
elseif z_use_s && z_use_u
    for k = 1:N-1
        grad_s_Zlin(:,k) = Cz(:,:,k)'*qz(:,k);
        grad_u_Zlin(:,k) = Dz(:,:,k)'*qz(:,k);
    end
elseif z_use_su && z_use_u
    for k = 1:N-1
        grad_su_Zlin(:,k) = Dsuz(:,:,k)'*qz(:,k);
        grad_u_Zlin(:,k) = Dz(:,:,k)'*qz(:,k);
    end
elseif z_use_s
    for k = 1:N-1
        grad_s_Zlin(:,k) = Cz(:,:,k)'*qz(:,k);
    end
elseif z_use_su
    for k = 1:N-1
        grad_su_Zlin(:,k) = Dsuz(:,:,k)'*qz(:,k);
    end
elseif z_use_u
    for k = 1:N-1
        grad_u_Zlin(:,k) = Dz(:,:,k)'*qz(:,k);
    end
end

if z_use_ter
    grad_s_Zlin_ter(:) = Cz_ter'*qz_ter;
end

end
