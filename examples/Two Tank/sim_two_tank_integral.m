clear yk rk tk h1k
r = [0.7;0];
h1 = x_prev(1);
h2 = x_prev(2);
e = 0;

tic
for k = 1:1000

if k>500
    r=[0.25;0];
end

y = C*x_prev ;


[u_prev,J,x0,t] = mpc_solve(x0,x_prev,u_prev,r,r(1),mpc,1e-2);
mpc.t = t;

e = e + Ts*(r(1)-h2);
h1 = h1 + Ts*(u_prev-sqrt(2*g)*sqrt(h1));
h2 = h2 + Ts*(sqrt(2*g)*sqrt(h1)-sqrt(2*g)*sqrt(h2));


x_prev = [h1;h2;e];


yk(1,k) = y(1);
rk(:,k) = r;
h1k(:,k) = h1;
tk(k) = t;
Jk(k) = J;
uk(k) = u_prev;

end
toc/1000
1/ans
J

close all
for y = 1:1
    figure
    plot(1:k,yk(y,:),1:k,rk(y,:),1:k,h1k(y,:),'g')
    grid on
end
legend('y(h2)','r','h1')
ylim([0 1])
figure
plot(1:k,tk)
grid on
figure
plot(1:k,uk)
hold on
plot(1:k-1,diff(uk))
grid on