clear yk rk tk h1k h2k
r = 0.7;
h1 = x_prev(1);
h2 = x_prev(2);

xf = h2;
tau = 0.1;

t1 = cputime;
for k = 1:500

if k>250
    r=0.25;
end

y = C*x_prev ;

xf = xf + Ts*(-xf/tau+r/tau);
xref = [xf;xf];

[u_prev,J,x0] = mpc_solve(x0,x_prev,u_prev,xf,[],mpc,xref,[],[]);

h1 = h1 + Ts*(u_prev-sqrt(2*g)*sqrt(h1));
h2 = h2 + Ts*(sqrt(2*g)*sqrt(h1)-sqrt(2*g)*sqrt(h2));

x_prev = [h1;h2];

yk(:,k) = xf;
rk(:,k) = r;
h1k(:,k) = h1;
h2k(:,k) = h2;
Jk(k) = J;
uk(k) = u_prev;

end
t2 = cputime-t1;
t2/500
500/t2

close all

figure
plot(Ts*(1:k),rk,'r',Ts*(1:k),h1k,'g',Ts*(1:k),h2k,'b')
grid on

legend('ref','h1','h2')
ylim([0 1])
figure
plot(Ts*(1:k),uk,Ts*(1:k-1),diff(uk))
grid on
legend('U','Delta U')