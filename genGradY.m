function [gradYmin,gradYmax] = genGradY(C,D,N,Nx,Nu,Ny,nx,nu,ny)
%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T

% y_min - C*x - D*sum_1_k(delta_u_k) - D*u0 <= 0 
% C*x + D*sum_1_k(delta_u_k) + D*u0 - y_max <= 0 

gradYmax = zeros(Nx+Nu,Ny);
for k = 1:N
    
    gradYmax(nu + nx*(k-1) + nu*(k-1) + 1 : nu + nx*(k)+ nu*(k),ny*(k-1)+1:ny*(k)) = [ C' ; D'];

end

gradYmin = -gradYmax;
end