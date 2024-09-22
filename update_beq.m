function beq = update_beq(beq,A,x_prev,N,nx,Bd,d)

for k = 1:N+1

    switch k 
        case 1
            beq((k-1)*nx+1:k*nx) = -A*x_prev-Bd*d(k);            
        otherwise

            beq((k-1)*nx+1:k*nx) = -Bd*d(k);
    end
end

end