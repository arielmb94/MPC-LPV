%% Cost function

e =R-Y;
Je = e'*(e.*Qe);

Ju = deltaU'*(deltaU.Qu);

JN_1 = x(Nx-nx+1:Nx)'*P*x(Nx-nx+1:Nx);

J = Je+Ju+JN_1;


%% Constraints

% deltaU
deltaU_Min-deltaU;
deltaU-deltaU_Max;

% U
U = zeros(Nu,1);
U(1:nu)=u0+deltaU(1:nu);
for i = 2:N

U(nu*(i-1)+1:nu*i) = U(nu*(i-2)+1:nu*(i-1))+deltaU(nu*(i-1)+1:nu*i);

end

U_min-U;
U-U_max;

% Y
Y = zeros(Ny,1);
for i = 1:N
i
    Y(ny*(i-1)+1:ny*i) = C*x(nx*(i-1)+1:nx*i)+D*U(nu*(i-1)+1:nu*i)+Dd*d(i+1);

end

Y_min-Y;
Y-Y_max;


% X

X_min-X(1:Nx-nx);
X(1:Nx-nx)-X_max;

% Terminal Set
X_min-X(1:Nx-nx);
X(1:Nx-nx)-X_max;