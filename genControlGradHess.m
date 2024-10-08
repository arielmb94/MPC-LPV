function [gradCtlrR,hessCtrlTerm] = genControlGradHess(R,N,Nx,Nu,nx,nu)


gradCtlrR = zeros(Nx+Nu,Nu);
gradCtlr = zeros(Nx+Nu,Nu);
%hessCtrlTerm = zeros(Nx+Nu,Nx+Nu);
for k = 0:N

    switch k
        case 0
            % Write in uo
            gradCtlrR(1:nu,1:nu) = R;
            gradCtlr(1:nu,1:nu) = eye(nu);
        otherwise
            % write in uk-1
            gradCtlrR(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k)+1:nu*(k+1)) = ...
                -R;
            gradCtlr(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k)+1:nu*(k+1)) = ...
                -eye(nu);
            % write in uk
            gradCtlrR(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k), nu*(k)+1:nu*(k+1)) = ...
                R;
            gradCtlr(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k), nu*(k)+1:nu*(k+1)) = ...
                eye(nu);

    end
end

%Grad term is then DeltaJu = gradCtlrR*deltaU

hessCtrlTerm = gradCtlrR*gradCtlr';


end