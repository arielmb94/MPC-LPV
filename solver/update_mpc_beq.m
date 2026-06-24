function mpc = update_mpc_beq(mpc,x_prev,d)

    if isfield(mpc, 'is_CT') && mpc.is_CT
        Ts = mpc.Ts;
        M_x0 = -eye(mpc.nx); % M=I
        Bd_term = Ts * mpc.Bd;
    else
        M_x0 = -mpc.A;
        Bd_term = mpc.Bd;
    end

    if isempty(d) || isempty(mpc.Bd)
        mpc.beq(1:mpc.nx) = M_x0 * x_prev;
    else
        for k = 0:mpc.N-1        
            switch k 
                case 0
                    mpc.beq(1:mpc.nx) = M_x0 * x_prev - Bd_term * d(1:mpc.nd);            
                otherwise        
                    mpc.beq(k*mpc.nx+1:(k+1)*mpc.nx) = ...
                                        -Bd_term * d(k*mpc.nd+1:(k+1)*mpc.nd);
            end
        end
    end

end