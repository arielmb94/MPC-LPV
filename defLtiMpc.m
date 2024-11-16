function mpc = defLtiMpc(N,A,B,C,D,Bd,Dd,Qe,Rdu,Ru,Cz,Dz,Ddz,Qz,Ci,Di,Ddi,...
    x_min,x_max,x_ter_min,x_ter_max,u_min,u_max,du_min,du_max,y_min,y_max,yi_min,yi_max)

mpc.N = N;              %prediction horizon

mpc.Qe = Qe;
mpc.Rdu = Rdu;
mpc.Ru = Ru;

mpc.A = A;
mpc.B = B;
mpc.Bd = Bd;
mpc.C = C;
mpc.D = D;
mpc.Dd = Dd;
% Performance Cost Matrix
mpc.Cz = Cz;
mpc.Dz = Dz;
mpc.Ddz = Ddz;
% General Inequality Matrix
mpc.Ci = Ci;
mpc.Di = Di;
mpc.Ddi = Ddi;

mpc.nx = size(A,1);  %number of states
mpc.nu = size(B,2);  %number of control inputs
mpc.nd = size(Bd,2);  %number of disturbance inputs
mpc.ny = size(C,1);  %number of measurements

mpc.ndz = size(Ddz,2);  %number of disturbance inputs to performance cost
mpc.nz = size(Cz,1);  %number of performances

mpc.ndyi = size(Ddi,2);  %number of disturbance inputs to general inequalities
mpc.nyi = size(Ci,1);  %number of general inequalities

mpc.Nx = (N+1)*mpc.nx;
mpc.Nu = (N+1)*mpc.nu;
mpc.Ny = (N+1)*mpc.ny;
mpc.Nd = (N+1)*mpc.nd;

mpc.Nz = N*mpc.nz;
mpc.Ndz = N*mpc.ndz;

mpc.Nyi = N*mpc.nyi;
mpc.Ndyi = N*mpc.ndyi;

mpc.x_min = x_min;
mpc.x_max = x_max;
mpc.x_ter_min = x_ter_min; 
mpc.x_ter_max = x_ter_max;
mpc.u_min = u_min;
mpc.u_max = u_max;
mpc.du_min = du_min;
mpc.du_max = du_max;
mpc.y_min = y_min;
mpc.y_max = y_max;
mpc.yi_min = yi_min;
mpc.yi_max = yi_max;

% A equality contraint (b equality constraints depends on x0 and d(k)
mpc.Aeq = genEqualities(A,B,N,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);
mpc.beq = zeros(size(mpc.Aeq,1),1);

% Cost Function Terms
mpc.hessCost = zeros(mpc.Nu+mpc.Nx);

if ~isempty(Qe)
    [mpc.gradErrQe,mpc.hessErrTerm] = genErrorGradHess(Qe,C,D,N,...
        mpc.Nx,mpc.Nu,mpc.Ny,mpc.nx,mpc.nu,mpc.ny);
    mpc.hessCost = mpc.hessCost + mpc.hessErrTerm;
else
    mpc.gradErrQe = [];
end

if ~isempty(Rdu)
    [mpc.gradDiffCtlrRdu,mpc.hessDiffCtrlTerm] = genDiffControlGradHess(Rdu,N,...
        mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);
    mpc.hessCost = mpc.hessCost + mpc.hessDiffCtrlTerm;
else
    mpc.gradDiffCtlrRdu = [];
end

if ~isempty(Ru)
    [mpc.gradCtlrRu,mpc.hessCtrlTerm] = genControlGradHess(Ru,N,...
        mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);
    mpc.hessCost = mpc.hessCost + mpc.hessCtrlTerm;
else
    mpc.gradCtlrRu = [];
end

if ~isempty(Qz)
    [mpc.gradPerfQz,mpc.hessPerfTerm] = genPerfGradHess(Qz,Cz,Dz,N,...
        mpc.Nx,mpc.Nu,mpc.Nz,mpc.nx,mpc.nu,mpc.nz);
    mpc.hessCost = mpc.hessCost + mpc.hessPerfTerm;
else
    mpc.gradPerfQz = [];
end

% Box Constraint Terms
m = 0; %constraint counter

% State box constraints
if ~isempty(mpc.x_min) || ~isempty(mpc.x_max)
    [mpc.gradXmin,mpc.gradXmax] = genGradX(N,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

    if ~isempty(mpc.x_min)
        [mpc.hessXmin,mi] = genHessIneq(mpc.gradXmin);
        m = m+mi;
    end
    if ~isempty(mpc.x_max)
        [mpc.hessXmax,mi] = genHessIneq(mpc.gradXmax);
        m = m+mi;
    end
end

% Terminal State box constraints
if ~isempty(mpc.x_ter_min) || ~isempty(mpc.x_ter_max)
    [mpc.gradXtermin,mpc.gradXtermax] = genGradXter(N,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

    if ~isempty(mpc.x_ter_min)
        [mpc.hessXtermin,mi] = genHessIneq(mpc.gradXtermin);
        m = m+mi;
    end
    if ~isempty(mpc.x_ter_max)
        [mpc.hessXtermax,mi] = genHessIneq(mpc.gradXtermax);
        m = m+mi;
    end
end

% Control box constraints
if ~isempty(mpc.u_min) || ~isempty(mpc.u_max)
    [mpc.gradUmin,mpc.gradUmax] = genGradU(N,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

    if ~isempty(mpc.u_min)
        [mpc.hessUmin,mi] = genHessIneq(mpc.gradUmin);
        m = m+mi;
    end
    if ~isempty(mpc.u_max)
        [mpc.hessUmax,mi] = genHessIneq(mpc.gradUmax);
        m = m+mi;
    end
end

% Differential Control box constraints
if ~isempty(mpc.du_min) || ~isempty(mpc.du_max)
    [mpc.gradDeltaUmin,mpc.gradDeltaUmax] = genGradDeltaU(N,mpc.Nx,mpc.Nu,mpc.nx,mpc.nu);

    if ~isempty(mpc.du_min)
        [mpc.hessDeltaUmin,mi] = genHessIneq(mpc.gradDeltaUmin);
        m = m+mi;
    end
    if ~isempty(mpc.du_max)
        [mpc.hessDeltaUmax,mi] = genHessIneq(mpc.gradDeltaUmax);
        m = m+mi;
    end
end

% Outputs box constraints
if ~isempty(mpc.y_min) || ~isempty(mpc.y_max)
    [mpc.gradYmin,mpc.gradYmax] = genGradY(C,D,N,mpc.Nx,mpc.Nu,mpc.Ny,...
        mpc.nx,mpc.nu,mpc.ny);

    if ~isempty(mpc.y_min)
        [mpc.hessYmin,mi] = genHessIneq(mpc.gradYmin);
        m = m+mi;
    end
    if ~isempty(mpc.y_max)
        [mpc.hessYmax,mi] = genHessIneq(mpc.gradYmax);
        m = m+mi;
    end
end

% General Inequalites box constraints
if ~isempty(mpc.yi_min) || ~isempty(mpc.yi_max)
    [mpc.gradYimin,mpc.gradYimax] = genGradY(Ci,Di,N,mpc.Nx,mpc.Nu,mpc.Nyi,...
        mpc.nx,mpc.nu,mpc.nyi);

    if ~isempty(mpc.yi_min)
        [mpc.hessYimin,mi] = genHessIneq(mpc.gradYimin);
        m = m+mi;
    end
    if ~isempty(mpc.yi_max)
        [mpc.hessYimax,mi] = genHessIneq(mpc.gradYimax);
        m = m+mi;
    end
end

mpc.m = m;


end