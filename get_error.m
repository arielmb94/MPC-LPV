function error = get_error(r,y,N,Ny)

    if length(r)<Ny

        r = repmat(r,N+1,1);

    end    

    error = r-y;
        
end