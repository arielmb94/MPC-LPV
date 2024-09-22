%% Generate Equality Constraint Matrix

%optimization variables are stored as [u0 x1 u1 x2 ... xn un xn+1]^T


Aeq = zeros(Nx,Nx+Nu);
beq = zeros(Nx,1);
for k = 1:N+1

    switch k 

        case 1
            Aeq((k-1)*nx+1:k*nx,1:nu+nx) = [B -eye(nx)];
            beq((k-1)*nx+1:k*nx) = -A*x_prev-Bd*d(k);
            
        otherwise

            Aeq((k-1)*nx+1:k*nx,(nu+nx)*(k-1)+1-nx:(nu+nx)*k) = [A B -eye(nx)];
            beq((k-1)*nx+1:k*nx) = -Bd*d(k);

    end


end
