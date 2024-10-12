function beq = update_beq(beq,A,x_prev,N,nx,Bd,d,Nd)

    if isempty(d) || isempty(Bd)

        beq(1:nx) = -A*x_prev;

    else

        % copy d N+1 times
        if length(d)<Nd
            d = repmat(d,N+1,1);
        end

        for k = 1:N+1
        
            switch k 
                case 1
                    beq((k-1)*nx+1:k*nx) = -A*x_prev-Bd*d(k);            
                otherwise
        
                    beq((k-1)*nx+1:k*nx) = -Bd*d(k);
            end
        end
    end

end