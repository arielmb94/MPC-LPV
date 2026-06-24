function mpc = genEqualities(mpc,A,B,N,N_h_ctr,nx,nu)

    % Backward euler
    if isfield(mpc, 'is_CT') && mpc.is_CT
        Ts = mpc.Ts;
        % Considering M = I (no mass matrix)
        M_s  = eye(nx);
        M_u  = Ts * B;
        M_sp = -(eye(nx) - Ts * A);
    else
        % Discrete time original system
        M_s  = A;
        M_u  = B;
        M_sp = -eye(nx);
    end

    for k = 0:N-1
        if k == 0
            mpc.Aeq(1:nx,1:nu+nx) = [M_u, M_sp];
        
        elseif k < N_h_ctr
            mpc.Aeq(k*nx+1:(k+1)*nx,(nu+nx)*k+1-nx:(nu+nx)*(k+1)) = [M_s, M_u, M_sp];

        else
            mpc.Aeq(k*nx+1:(k+1)*nx,(nx+nu)*(N_h_ctr-1)+1:(nx+nu)*(N_h_ctr-1)+nu) = ...
                M_u;

            mpc.Aeq(k*nx+1:(k+1)*nx,(nx+nu)*(N_h_ctr-1)+nu+nx*(k-(N_h_ctr))+1:...
                (nx+nu)*(N_h_ctr-1)+nu+nx*(k-(N_h_ctr-2))) = ...
                [M_s, M_sp];

        end
    end
end