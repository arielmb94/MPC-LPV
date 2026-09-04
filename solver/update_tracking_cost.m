function [update_tracking,recompute_cost_hess,...
          grad_u_E_0,grad_s_E,grad_u_E,grad_s_E_ter,...
          R_E_0,Q_E,R_E,Y_E,Q_E_ter] = update_tracking_cost(...
          update_tracking,recompute_cost_hess,...
          grad_u_E_0,grad_s_E,grad_u_E,grad_s_E_ter,...
          R_E_0,Q_E,R_E,Y_E,Q_E_ter,...
          y_use_k0,y_use_s,y_use_u,y_use_ter,...
          Qe_0,Qe,Qe_ter,D_0,C,D,C_ter,N)

update_tracking = 0;
recompute_cost_hess = 1;

if y_use_k0
    grad_u_E_0(:,:) = D_0'*Qe_0;
    R_E_0(:,:) = grad_u_E_0*D_0;
end

if y_use_s && y_use_u
    for k = 1:N-1
        grad_s_E(:,:,k) = C(:,:,k)'*Qe(:,:,k);
        grad_u_E(:,:,k) = D(:,:,k)'*Qe(:,:,k);

        Q_E(:,:,k) = grad_s_E(:,:,k)*C(:,:,k);
        R_E(:,:,k) = grad_u_E(:,:,k)*D(:,:,k);
        Y_E(:,:,k) = grad_u_E(:,:,k)*C(:,:,k);
    end
elseif y_use_s
    for k = 1:N-1
        grad_s_E(:,:,k) = C(:,:,k)'*Qe(:,:,k);
        Q_E(:,:,k) = grad_s_E(:,:,k)*C(:,:,k);
    end
elseif y_use_u
    for k = 1:N-1
        grad_u_E(:,:,k) = D(:,:,k)'*Qe(:,:,k);
        R_E(:,:,k) = grad_u_E(:,:,k)*D(:,:,k);
    end
end

if y_use_ter
    grad_s_E_ter(:,:) = C_ter'*Qe_ter;
    Q_E_ter(:,:) = grad_s_E_ter*C_ter;
end

end
