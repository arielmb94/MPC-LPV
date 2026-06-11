function mpc = update_A_Equalities(mpc)

% k = 0
%Dynamics [B -I]*[u0 s1]'=-A*s0-Dd*d0

for k = 2:mpc.N
    %Dynamics [A B -I]*[s1 u1 s2]'
    row = mpc.dyn_k(:,k);
    s_col = mpc.s_index_k(:,k);
    mpc.Aeq(row,s_col) = mpc.A;
end
    
end



