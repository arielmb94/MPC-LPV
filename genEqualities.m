function Aeq = genEqualities(A,B,N,Nx,Nu,nx,nu)

Aeq = zeros(Nx,Nx+Nu);
for k = 1:N+1
    switch k 
        case 1
            Aeq((k-1)*nx+1:k*nx,1:nu+nx) = [B -eye(nx)];
        otherwise
            Aeq((k-1)*nx+1:k*nx,(nu+nx)*(k-1)+1-nx:(nu+nx)*k) = [A B -eye(nx)];
    end
end

end