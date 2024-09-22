%% Example random system
sys = drss(4,3,3)

A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

Bd = B(:,end);
B = B(:,1:end-1);

Dd = D(:,end);
D = D(:,1:end-1);

%%
% Define mpc structure



