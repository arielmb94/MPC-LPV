clear yk rk tk h1k
r = 0.7*ones(ny,1);
h1 = x_prev(1);
h2 = x_prev(2);

tic
for k = 1:500

if k>250
    r=0.45*ones(ny,1);
end

y = C*x_prev ;


[u_prev,J,x0,t] = mpc_solve(x0,x_prev,u_prev,r,d,mpc,1e-2,1);
mpc.t = t;

h1 = h1 + Ts*(u_prev-sqrt(2*g)*sqrt(h1));
h2 = h2 + Ts*(sqrt(2*g)*sqrt(h1)-sqrt(2*g)*sqrt(h2));

x_prev = [h1;h2];


yk(:,k) = y;
rk(:,k) = r;
h1k(:,k) = h1;
tk(k) = t;
Jk(k) = J;

end
toc/500
1/ans
J

close all
for y = 1:ny
    figure
    plot(1:k,yk(y,:),1:k,rk(y,:),1:k,h1k(y,:),'g')
    grid on
end
legend('y(h2)','r','h1')
ylim([0 1])
figure
plot(1:k,tk)
grid on