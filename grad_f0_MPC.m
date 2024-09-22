function grad_J = grad_f0_MPC(Nx,Nu,gradErrQe,err,gradCtlrR,deltaU,P,s_ter,nx)

    grad_J = zeros(Nx+Nu,1);

    if ~isempty(gradErrQe)

        grad_J = grad_J - gradErrQe*err;

    end    

    if ~isempty(gradCtlrR)

        grad_J = grad_J + gradCtlrR*deltaU;

    end  

    if ~isempty(P)

        gradTer = P*s_ter;

        grad_J(end-nx+1:end) = grad_J(end-nx+1:end) + gradTer;

    end  

        
end



