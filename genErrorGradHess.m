function [gradErrQe,hessErrTerm] = genErrorGradHess(Qe,C,D,N,Nx,Nu,Ny,nx,nu,ny)

%Gradient = gradError*Qe*error
gradErrQe = zeros(Nx+Nu,Ny);
hessErrTerm = zeros(Nx+Nu,Nx+Nu);
for k = 0:N-1

    switch k

        case 0
            % gradError = -D'*Qe
            % Grad term is then DeltaJe = -gradErrQe*err
            gradErrQe(1:nu, 1:ny) = D'*Qe;           
            % hessError = gradErrQe*D
            hessErrTerm(1:nu, 1:nu) = gradErrQe(1:nu, 1:ny)*D;

        otherwise
            % gradError = -[C';D']*Qe
            % Grad term is then DeltaJe = -gradErrQe*err
            gradErrQe(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),ny*(k)+1:ny*(k+1)) = ...
                [ C' ; D']*Qe;
            % hessError = gradErrQe*[C D]
            hessErrTerm(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k))...
                = gradErrQe(nu + 1 + nx*(k-1) + nu*(k-1)  : nu + nx*(k)+ nu*(k),ny*(k)+1:ny*(k+1))*[C D];

    end
    
end

end