function beq = update_beq(beq,A,x_prev,N,nx,Bd,d,Nd)

    if isempty(d) || isempty(Bd)

        beq(1:nx) = -A*x_prev;

    else

        % copy d N times
        if length(d)<Nd
            d = repmat(d,N,1);
        end

        for k = 0:N-1
        
            switch k 
                case 0
                    beq(1:nx) = -A*x_prev-Bd*d(k);            
                otherwise
        
                    beq(k*nx+1:(k+1)*nx) = -Bd*d(k);
            end
        end
    end

end