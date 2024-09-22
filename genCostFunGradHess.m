
iter_num = 1;
%% Error Term

%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

for i = 1:iter_num
tic
%Gradient = gradError*Qe*error
gradErrQe = zeros(Nx+Nu,Ny);
hessErrTerm = zeros(Nx+Nu,Nx+Nu);
for k = 1:N
    

    gradErrQe(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),ny*(k-1)+1:ny*(k)) = ...
        [ C' ; D']*Qe;
    hessErrTerm(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k))...
        = gradErrQe(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),ny*(k-1)+1:ny*(k))*[C D];

end
timing(i) = toc;

end
%Grad term is then DeltaJe = -gradErrQe*err




%% Delta Control Term

%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

for i = 1:iter_num
tic
gradCtlrR = zeros(Nx+Nu,Ny);
hessCtrlTerm = zeros(Nx+Nu,Nx+Nu);
for k = 0:N

    switch k
        case 0
            % Write in uo
            gradCtlrR(1:nu,1:nu) = R;
        otherwise
            % write in uk-1
            gradCtlrR(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k)+1:nu*(k+1)) = ...
                -R;
            % write in uk
            gradCtlrR(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k), nu*(k)+1:nu*(k+1)) = ...
                R;

    end
end

%Grad term is then DeltaJu = gradCtlrR*deltaU

hessCtrlTerm = gradCtlrR*gradDeltaUmax';

timing(i) = toc;
end



%% Terminal Cost

%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

gradXn1term = zeros(Nx+Nu,1);
%gradXn1term(end-nx+1:end) = P*x(end-nx+1:end);

HessXn1term = zeros(Nx+Nu,Nx+Nu);
HessXn1term(Nx+Nu-nx+1:end,Nx+Nu-nx+1:end) = P;


