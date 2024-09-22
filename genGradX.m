function [gradXmin,gradXmax] = genGradX(N,Nx,Nu,nx,nu)
%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% x_min - xk <= 0 
% xk - xmax <= 0 

gradXmax = zeros(Nx+Nu,Nx);
for k = 1:N+1
    gradXmax(nu + nx*(k-1)+ + nu*(k-1) + 1 : nu + nx*(k)+ + nu*(k-1), nx*(k-1)+1:nx*(k)) = eye(nx);
end

gradXmin = -gradXmax;
end