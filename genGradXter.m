function [gradXtermin,gradXtermax] = genGradXter(N,Nx,Nu,nx,nu)
%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% x_min - xk <= 0 
% xk - xmax <= 0 

gradXtermax = zeros(Nx+Nu,nx);
k = N;

gradXtermax(nu + nx*(k-1)+ nu*(k-1) + 1 : nu + nx*(k)+ nu*(k-1), :) = eye(nx);
gradXtermin = -gradXtermax;

end