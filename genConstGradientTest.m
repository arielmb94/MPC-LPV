%% Gradient for constraints on DeltaU

%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% u_min - uk + uk-1 <= 0 
% uk - uk-1 - umax <= 0 
tic

gradDeltaUmax = zeros(Nx+Nu,Nu);

for k = 0:N

    switch k
        case 0
            % Write in uo
            gradDeltaUmax(1:nu,1:nu) = eye(nu);
        otherwise
            % write in uk-1
            gradDeltaUmax(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k)+1:nu*(k+1)) = -eye(nu);
            % write in uk
            gradDeltaUmax(nu + nx*k + nu*(k-1)+1: nu + nx*k + nu*(k), nu*(k)+1:nu*(k+1)) = eye(nu);

    end
end

gradDeltaUmin = -gradDeltaUmax;

%% Gradient for constraints on U

%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% u_min - uk <= 0 
% uk - umax <= 0 

gradUmax = zeros(Nx+Nu,Nu);

for k = 1:N+1


    gradUmax(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k-1)+1:nu*(k)) = eye(nu);

end

gradUmin = -gradUmax;

%% Gradient for constraints on X 

%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% x_min - xk <= 0 
% xk - xmax <= 0 

gradXmax = zeros(Nx+Nu,Nx);
for k = 1:N+1

    gradXmax(nu + nx*(k-1)+ + nu*(k-1) + 1 : nu + nx*(k)+ + nu*(k-1), nx*(k-1)+1:nx*(k)) = eye(nx);

end

gradXmin = -gradXmax;
%% Gradient on terminal state constraints on X := xn+

%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

gradXn_1max = zeros(Nx+Nu,nx);
k = N+1;
gradXn_1max(nu + nx*(k-1)+ + nu*(k-1) + 1 : nu + nx*(k)+ + nu*(k-1), 1:nx) = eye(nx);

gradXn_1min = -gradXn_1max;


%% Gradients for constraints on Y

% y_min - C*x - D*sum_1_k(delta_u_k) - D*u0 <= 0 
% C*x + D*sum_1_k(delta_u_k) + D*u0 - y_max <= 0 

gradYmax = zeros(Nx+Nu,Ny);
for k = 1:N
    
    gradYmax(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),ny*(k-1)+1:ny*(k)) = [ C' ; D'];

end

gradYmin = -gradYmax;




