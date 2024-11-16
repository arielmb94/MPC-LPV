function [gradPerfQz,hessPerfTerm] = genPerfGradHess(Qz,C,D,N,Nx,Nu,Nz,nx,nu,nz)

%Gradient = gradPerf*Qz*perf
gradPerfQz = zeros(Nx+Nu,Nz);
hessPerfTerm = zeros(Nx+Nu,Nx+Nu);
for k = 1:N
    gradPerfQz(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),nz*(k-1)+1:nz*(k)) = ...
        [ C' ; D']*Qz;
    %Grad term is then DeltaJz = gradPerfQz*z

    hessPerfTerm(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k))...
        = gradPerfQz(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),nz*(k-1)+1:nz*(k))*[C D];
end

end