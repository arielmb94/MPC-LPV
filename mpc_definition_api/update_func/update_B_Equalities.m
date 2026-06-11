function mpc = update_B_Equalities(mpc)

% k = 0
%Dynamics [B -I]*[u0 s1]'=-A*s0-Dd*d0
row = mpc.dyn_k(:,1);
u_col = mpc.u_index_k(:,1);

mpc.Aeq(row,u_col) = mpc.B;

for k = 2:mpc.N
    %Dynamics [A B -I]*[s1 u1 s2]'
    row = mpc.dyn_k(:,k);
    u_col = mpc.u_index_k(:,k);

    mpc.Aeq(row,u_col) = mpc.B;
end
    
end



