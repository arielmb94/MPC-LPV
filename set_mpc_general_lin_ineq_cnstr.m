function mpc = set_mpc_general_lin_ineq_cnstr(mpc,yi_min,yi_max,Ci,Di,Ddi)

mpc.yi_min = yi_min;
mpc.yi_max = yi_max;

% General Inequality Matrix
mpc.Ci = Ci;
mpc.Di = Di;
mpc.Ddi = Ddi;

mpc.ndi = size(Ddi,2);  %number of disturbance inputs to general inequalities
mpc.nyi = size(Ci,1);  %number of general inequalities

if mpc.Di == 0
    mpc.Nyi = (N-1)*mpc.nyi;
else
    mpc.Nyi = N*mpc.nyi;
end

mpc.Ndi = N*mpc.ndi;

% General Inequalites box constraints
[mpc.gradYimin,mpc.gradYimax] = genGradY(mpc.Ci,mpc.Di,mpc.N,mpc.Nx,mpc.Nu,mpc.Nyi,...
    mpc.nx,mpc.nu,mpc.nyi);

if ~isempty(mpc.yi_min)
    [mpc.hessYimin,mi] = genHessIneq(mpc.gradYimin);
    mpc.m = mpc.m+mi;
end
if ~isempty(mpc.yi_max)
    [mpc.hessYimax,mi] = genHessIneq(mpc.gradYimax);
    mpc.m = mpc.m+mi;
end

end