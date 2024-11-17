function [s,s_all,s_ter] = get_x(x,s_prev,nx,nu,N,Nx)

    s = zeros(Nx,1);
    for k = 1:N+1
        
        s((k-1)*nx+1:k*nx) = x(nx*(k-1) + nu*k + 1 : nx*k + nu*k);
    
    end

    s_all = [s_prev;s];
    s_ter = x(nx*(N) + nu*(N+1) + 1 : end);
        
end