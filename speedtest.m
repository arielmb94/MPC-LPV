iter_num = 1e3;
clear timing
% test
%%
%timing = zeros(iter_num,2);
tic
for k = 1:iter_num

    mpc = defLtiMpc(N,A,B,C,D,Bd,Dd,Qe,R,x_min,x_max,u_min,u_max,du_min,du_max,y_min,y_max);


end
timing(1) = toc/iter_num
