function mpc = genEqualities(mpc)

du = mpc.has_du;
N = mpc.N

r_len = (mpc.Nx) + sum(mpc.ng_k);
c_len = mpc.Nx+mpc.Nu+mpc.Nu*du+sum(mpc.ng_k)+sum(mpc.nv_k);

% crear indice de cada inequality

Aeq = zeros(r_len,c_len);
beq = zeros(r_len,1);

% k = 0
%Dynamics [B 0 -I 0]*[u0 g0 s1 su1]'
row = mpc.dyn_k(:,1);
col = mpc.u_index_k(:,1);
Aeq(row,col) = mpc.B;
col = mpc.s_index_k(:,2);
Aeq(row,col) = -eye(mpc.nx);
%[I 0 0 -I]*[u0 g0 s1 su1]'
% row = mpc.nx+1:mpc.nu*du;
% col = mpc.u_index_k(:,1);
% Aeq(row,col) = eye(nu*du);
% col = mpc.u_index_k(:,1);
% Aeq(row,col) = -eye(nu*du);

if mpc.has_u_cnstr
    if mpc.u_cnstr.min_limit
        %Ineq [-I I]*[u0 g0]'
        row = mpc.u_cnstr.min_eq_index_k(:,1);
        col = mpc.u_index_k(:,1);
        Aeq(row,col) = -eye(mpc.nu);
        col = mpc.u_cnstr.g_min_index_k(:,1);
        Aeq(row,col) = eye(mpc.nu);

        beq(row) = -mpc.u_cnstr.min;
    end
    if mpc.u_cnstr.max_limit
        %Ineq [I I]*[u0 g0]'
        row = mpc.u_cnstr.max_eq_index_k(:,1);
        col = mpc.u_index_k(:,1);
        Aeq(row,col) = eye(mpc.nu);
        col = mpc.u_cnstr.g_max_index_k(:,1);
        Aeq(row,col) = eye(mpc.nu);

        beq(row) = mpc.u_cnstr.max;
    end
end


for k = 2:N
    %Dynamics [A B -I]*[s1 u1 s2]'
    row = mpc.dyn_k(:,k);
    col = mpc.s_index_k(:,k);
    Aeq(row,col) = mpc.A;
    col = mpc.u_index_k(:,k);
    Aeq(row,col) = mpc.B;
    col = mpc.s_index_k(:,k+1);
    Aeq(row,col) = -eye(mpc.nx);

    if mpc.has_s_cnstr
        if mpc.s_cnstr.min_limit
            %Ineq [-I I -I]*[sk gk vk]'
            row = mpc.s_cnstr.min_eq_index_k(:,k);
            col = mpc.s_index_k(:,k);
            Aeq(row,col) = -eye(mpc.nx);
            col = mpc.s_cnstr.g_min_index_k(:,k);
            Aeq(row,col) = eye(mpc.nx);
            col = mpc.s_cnstr.v_min_index_k(:,k);
            Aeq(row,col) = -eye(mpc.nx);

            beq(row) = -mpc.s_cnstr.min;
        end
        if mpc.s_cnstr.max_limit
            %Ineq [I I -I]*[sk gk vk]'
            row = mpc.s_cnstr.max_eq_index_k(:,k);
            col = mpc.s_index_k(:,k);
            Aeq(row,col) = eye(mpc.nx);
            col = mpc.s_cnstr.g_max_index_k(:,k);
            Aeq(row,col) = eye(mpc.nx);
            col = mpc.s_cnstr.v_max_index_k(:,k);
            Aeq(row,col) = -eye(mpc.nx);

             beq(row) = mpc.s_cnstr.max;
        end
    end

    if mpc.has_u_cnstr
        if mpc.u_cnstr.min_limit
            %Ineq [-I I]*[u0 g0]'
            row = mpc.u_cnstr.min_eq_index_k(:,k);
            col = mpc.u_index_k(:,k);
            Aeq(row,col) = -eye(mpc.nu);
            col = mpc.u_cnstr.g_min_index_k(:,k);
            Aeq(row,col) = eye(mpc.nu);

            beq(row) = -mpc.u_cnstr.min;
        end
        if mpc.u_cnstr.max_limit
            %Ineq [I I]*[u0 g0]'
            row = mpc.u_cnstr.max_eq_index_k(:,k);
            col = mpc.u_index_k(:,k);
            Aeq(row,col) = eye(mpc.nu);
            col = mpc.u_cnstr.g_max_index_k(:,k);
            Aeq(row,col) = eye(mpc.nu);

            beq(row) = mpc.u_cnstr.max;
        end
    end

end

if mpc.has_s_cnstr
    if mpc.s_cnstr.min_limit
        %Ineq [-I I -I]*[sk gk vk]'
        row = mpc.s_cnstr.min_eq_index_k(:,N+1);
        col = mpc.s_index_k(:,N+1);
        Aeq(row,col) = -eye(mpc.nx);
        col = mpc.s_cnstr.g_min_index_k(:,N+1);
        Aeq(row,col) = eye(mpc.nx);
        col = mpc.s_cnstr.v_min_index_k(:,N+1);
        Aeq(row,col) = -eye(mpc.nx);

        beq(row) = -mpc.s_cnstr.min;
    end
    if mpc.s_cnstr.max_limit
        %Ineq [I I -I]*[sk gk vk]'
        row = mpc.s_cnstr.max_eq_index_k(:,N+1);
        col = mpc.s_index_k(:,N+1);
        Aeq(row,col) = eye(mpc.nx);
        col = mpc.s_cnstr.g_max_index_k(:,N+1);
        Aeq(row,col) = eye(mpc.nx);
        col = mpc.s_cnstr.v_max_index_k(:,N+1);
        Aeq(row,col) = -eye(mpc.nx);

        beq(row) = mpc.s_cnstr.max;
    end
end

mpc.Aeq = Aeq;
mpc.beq = beq;
    
end


function [A_k b_k] = genInequality(mpc,nu,nv,C,D,Duprev,Dd)
    
    A_k = [C Duprev D zeros(nm,nx) zeros(nm,nu) zeros(nm,nu) zeros(nm,nu) eye(nm,nh) ]
    b_k = -Dd;

end


function [A_k b_k] = genEqualities_k0(mpc,B)
    du = mpc.has_du;
    nx = mpc.nx;
    nu = mpc.nu;
    ng = mpc.ng_k(1);

    A_k_dyn = [B zeros(nx,ng) -eye(nx) zeros(nx,du*nu);
               eye(du*nu) zeros(du*nu,ng) zeros(nu,du*nx) -eye(du*nu)];

    if mpc.has_u_cnstr


    end
end