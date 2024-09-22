function [gradUmin,gradUmax] = genGradU(N,Nx,Nu,nx,nu)
%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% u_min - uk <= 0 
% uk - umax <= 0 

gradUmax = zeros(Nx+Nu,Nu);

for k = 1:N+1
    gradUmax(nu + nx*(k-1) + nu*(k-2)+1: nu + nx*(k-1) + nu*(k-1), nu*(k-1)+1:nu*(k)) = eye(nu);
end

gradUmin = -gradUmax;
end