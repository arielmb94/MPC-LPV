function [gradErrQe,hessErrTerm] = genErrorGradHess(Qe,C,D,N,Nx,Nu,Ny,nx,nu,ny)

%Gradient = gradError*Qe*error
gradErrQe = zeros(Nx+Nu,Ny);
hessErrTerm = zeros(Nx+Nu,Nx+Nu);
for k = 1:N
    gradErrQe(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),ny*(k-1)+1:ny*(k)) = ...
        [ C' ; D']*Qe;
    %Grad term is then DeltaJe = -gradErrQe*err

    hessErrTerm(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k))...
        = gradErrQe(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),ny*(k-1)+1:ny*(k))*[C D];
end

end