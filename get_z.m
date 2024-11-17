function z = get_z(s,u,d,nx,nu,nz,ndz,N,Nz,C,D,Dd,Ndz)

    if isempty(d) || isempty(Dd)

        z = zeros(Nz,1);
        for k = 1:N
    
            sk = s((k-1)*nx+1:k*nx);
            uk = u((k-1)*nu+1:k*nu);
    
            z((k-1)*nz+1:k*nz) = C*sk+D*uk;
        
        end
    
    else

        % copy d N+1 times
        if length(d)<Ndz
            d = repmat(d,N+1,1);
        end

        z = zeros(Nz,1);
        for k = 1:N
    
            sk = s((k-1)*nx+1:k*nx);
            uk = u((k-1)*nu+1:k*nu);
            dk = d((k-1)*ndz+1:k*ndz);
    
            z((k-1)*nz+1:k*nz) = C*sk+D*uk+Dd*dk;
        
        end
        
    end
end