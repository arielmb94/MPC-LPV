function mpc = genEqualities(mpc)

nx = mpc.nx;
nu = mpc.nu;
nse = mpc.nse;

s_col = 1:nx;
su_col = nx+1:nse;
u_col = nse+1:nse+nu;

mpc.s_col = s_col;
mpc.su_col = su_col;
mpc.se_col = [s_col su_col];
mpc.u_col = u_col;

mpc.rp_0 = zeros(nse,1);
mpc.rp_k = zeros(nse,mpc.N-1);
mpc.beq_0 = zeros(nx,1);
mpc.beq_k = zeros(nx,mpc.N-1);


%% dynamics

if mpc.has_du
    mpc.A_kkt = zeros(mpc.nse,mpc.nse,mpc.N-1);
    mpc.A_kkt(mpc.s_col,mpc.s_col,:) = mpc.A(:,:,2:mpc.N);

    mpc.B_kkt_0 = [mpc.B(:,:,1);eye(mpc.nu)];
    mpc.B_kkt = zeros(mpc.nse,mpc.nu,mpc.N-1);
    mpc.B_kkt(mpc.s_col,:,:) = mpc.B(:,:,2:mpc.N);
    mpc.B_kkt(mpc.su_col,:,:) = eye(mpc.nu).*ones(1,1,mpc.N-1);

else
    mpc.A_kkt = mpc.A(:,:,2:mpc.N);

    mpc.B_kkt_0 = mpc.B(:,:,1);
    mpc.B_kkt = mpc.B(:,:,2:mpc.N);
end

end
