clear yk rk
r = -2*ones(ny,1);
tic
for k = 1:100

dist = randn*0.1;

if k>50
    r=2*ones(ny,1);
end

y = C*x_prev + D*u_prev + Dd*dist;


[u_prev,J,x0] = mpc_solve(x0,x_prev,u_prev,r,d,mpc,1e-2,0.01);

x_prev = A*x_prev + B*u_prev + Bd*dist;


yk(:,k) = y;
rk(:,k) = r;
Jk(k) = J;

end
toc/100
1/ans
J

close all
for y = 1:ny
    figure
    plot(1:k,yk(y,:),1:k,rk(y,:))
    grid on
end


%%

%states
s = get_x(x0,mpc.nx,mpc.nu,mpc.N,mpc.Nx)
% control actions
u = get_u(x0,mpc.nx,mpc.nu,mpc.N,mpc.Nu)
u0 = u(1:mpc.nu)
% differential control action
du = diff_u(u,u_prev,mpc.nu,mpc.N,mpc.Nu)
% system outputs
y = get_y(s,u,d,mpc.nx,mpc.nu,mpc.ny,mpc.nd,mpc.N,mpc.Ny,...
    mpc.C,mpc.D,mpc.Dd)
% error signal
err = get_error(r,y,mpc.N,mpc.Ny)





