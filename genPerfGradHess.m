function [gradPerfQz,hessPerfTerm] = genPerfGradHess(Qz,C,D,N,Nx,Nu,Nz,nx,nu,nz)

%Gradient = gradPerf*Qz*perf
gradPerfQz = zeros(Nx+Nu,Nz);
hessPerfTerm = zeros(Nx+Nu,Nx+Nu);
for k = 0:N-1

    switch k

        case 0
            gradPerfQz(1 : nu, 1:nz) =  D'*Qz;
            %Grad term is then DeltaJz = gradPerfQz*z
        
            hessPerfTerm(1 : nu , 1 : nu) = gradPerfQz(1 : nu)*D;

        otherwise
            gradPerfQz(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),nz*(k)+1:nz*(k+1)) = ...
                [ C' ; D']*Qz;
            %Grad term is then DeltaJz = gradPerfQz*z
        
            hessPerfTerm(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k))...
                = gradPerfQz(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),nz*(k)+1:nz*(k+1))*[C D];
    end
end

end