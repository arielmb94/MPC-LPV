function y = get_y(s,u,d,nx,nu,ny,nd,N,Ny,C,D,Dd,Nd)

    if isempty(d) || isempty(Dd)

        y = zeros(Ny,1);
        for k = 1:N
    
            sk = s((k-1)*nx+1:k*nx);
            uk = u((k-1)*nu+1:k*nu);
    
            y((k-1)*ny+1:k*ny) = C*sk+D*uk;
        
        end
    
    else

        % copy d N+1 times
        if length(d)<Nd
            d = repmat(d,N+1,1);
        end

        y = zeros(Ny,1);
        for k = 1:N
    
            sk = s((k-1)*nx+1:k*nx);
            uk = u((k-1)*nu+1:k*nu);
            dk = d((k-1)*nd+1:k*nd);
    
            y((k-1)*ny+1:k*ny) = C*sk+D*uk+Dd*dk;
        
        end
        
    end
end