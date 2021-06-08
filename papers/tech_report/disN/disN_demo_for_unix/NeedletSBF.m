function NeedletSBF(J)
% Test evaluation of Needlet Psi function

format compact

% subdir = '\'; % for windows
subdir = '/'; % for linux

sv_dir = ['disN_rbf' subdir];
if ~exist(sv_dir,'dir')
    mkdir(sv_dir);
end
fprintf('\n******************************************\n');
% dimension of the sphere
d = 2;
fprintf('Block algorithm for discrete needlet approxiamtions for RBFs on S^%d', d);

% Symmetric spherical designs for needlet and discretization quadratures
% ip = 17;
QH = 'SD'; %%%%%%
QN = 'SD';
QNm = 'SD';
funtxt_ce = {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'};

% Maximum needlet level
% J = 6; %%%%%%

% Area of sphere for noamalized surface measure
% S2area = 4*pi;
S2area = 1;

% Quadrature for for discretizing inner product exact for degree degJ
degJ = 3*2^(J-1) - 1;
[W_QH, Y_QH] = pntset(QH, degJ);
W_QH = W_QH/S2area;
NY = size(Y_QH,2);
% Evaluate needlet approximation at points in Y to estimate errors
deg_Nm = 301; %%%%%%
[~,x_Nm] = pntset(QNm, deg_Nm);
N = size(x_Nm,2);
fprintf('\n******************************************\n');
fprintf('Discretisation quadrature rule: QH = %s', QH);
fprintf(', J = %d, deg = %d, NY = %d\n', J, degJ, NY);
fprintf('Evaluating points: QNm = %s', QNm);
fprintf(', N = %d\n', N);
% Filter
kappa = 5;
h = @(t) Needletfilterpoly(kappa, t);
fprintf('Needlet filter with smoothnss kappa = %d', kappa);
fprintf('\n******************************************\n');
fprintf('Test functions: normalised Wendland f_k, ');
fprintf('k = 0,1,2,3,4');
fprintf('\n******************************************\n');

for i_ce = 1:length(funtxt_ce)
    funtxt = funtxt_ce{i_ce};
     switch funtxt
         case 'nrWend_k0'
             rbf_k = 0;             
         case 'nrWend_k1'
             rbf_k = 1;
         case 'nrWend_k2'
             rbf_k = 2;
         case 'nrWend_k3'
             rbf_k = 3;
         case 'nrWend_k4'
             rbf_k = 4;
     end
     switch funtxt
         case {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'}
             fun = @rbf_nr;
     end

fprintf('\n====== Approximating for f_{%d} starts ======\n',rbf_k);
fprintf('Discretisation quadrature rule: QH = %s', QH);
fprintf(', J = %d, deg = %d, NY = %d\n', J, degJ, NY);

% Test function at discretization points
f_QH = fun(Y_QH',rbf_k)';
switch funtxt
    case {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'}
        fun_fprt = 'normalised Wendland';
        k = rbf_k;
        delta = (3*k+3)*gamma(k+1/2)/(2*gamma(k+1));
end
fprintf('Test function: %s, k = %d, delta = %.2f\n', fun_fprt, k, delta);
% Weighted function values for discretized inner product
f_QH = f_QH.*W_QH;
f_QH = f_QH(:);
% clear W

qN = cell(J+1,1);
t00 = tic;
% Evaluate discrete inner product for needlet coefficients at all levels
for j = 0:J
    % Needlet quadrature rule exact up to degree 2^(j+1)-1
    deg = 2^(j+1)-1;
    [w_QN, X_QN] = pntset(QN, deg);
    Nj = size(X_QN,2);
    w_QN = w_QN/S2area;
    
    % Block algorithm to evaluate all needlets of level j at points in Y
    t0 = tic;
    blksz = 1020;
    nblk = ceil(NY/blksz);
    PsiY = zeros(Nj, NY);
    for  blk = 1:nblk
        % Calculate Needlet values for points Y(:,I)
        I = (blk-1)*blksz+1: min(blk*blksz,NY);
        PsiY(:,I) =  NeedletPsi(j, X_QN, w_QN, h, Y_QH(:,I));
    end
    t_psi = toc(t0);
    fprintf('Evaluating discrete inner products with needlets psi_{jk}, k = 1,...,Nj at N points\n');
    fprintf('Order j = %d, Nj = %d, N = %d, ', j, Nj, NY);
    fprintf('BlkSz = %d, time = %.2f secs\n', blksz, t_psi);    
    qN{j+1} = PsiY*f_QH;    
end
clear PsiY
tQN = toc(t00);
fprintf(' == Evalauting discrete needlet coefficients for J = %d, N = %d', J, NY);
fprintf(', Time = %.2f secs\n', tQN);

% Evaluate needlet approximation at points in Y to estimate errors
% deg_Nm = 301; %%%%%%
% [~,x_Nm] = pntset(QNm, deg_Nm);
% N = size(x_Nm,2);
fprintf('Evaluating points: QNm = %s', QNm);
fprintf(', N = %d\n', N);

% Evaluate approxiamtion all all points in Y
t00 = tic;
U = zeros(1,N);
dN = zeros(J+1,N);
for j = 0:J
    
    % Needlet quadrature rule exact up to degree 2^(j+1)-1
    deg = 2^(j+1)-1;
    %deg = 2^(j+1)+1;
    [w_QN, X_QN] = pntset(QN, deg);
    Nj = size(X_QN,2);
    w_QN = w_QN/S2area;
    
    % Needlet coefficients for level j
    QNj = qN{j+1}';
    
    % Block algorithm to evaluate all needlets of level j at points in Y
    t0 = tic;
    % Choose blksz to balance storage and speed
%     blksz = 1020;
    blksz = max(8e2,floor(N/Nj))+1;
    nblk = ceil(N/blksz);
    % Evaluating all needlets of level j for all points is more efficient
    % but requires a lot of storage
    for  blk = 1:nblk
        % Calculate Needlet values for points Y(:,I)
        I = (blk-1)*blksz+1: min(blk*blksz,N);
        PsiYI =  NeedletPsi(j, X_QN, w_QN, h, x_Nm(:,I));
        U(I) = U(I) + QNj*PsiYI;
    end
    
    dN(j+1,:) = U;
    
    t_psi = toc(t0);
    fprintf('Evaluating needlet approximation with psi_{jk}, k = 1,...,Nj at N points\n');
    fprintf('Order j = %d, Nj = %d, N = %d, ', j, Nj, N);
    fprintf('BlkSz = %d, time = %.2f secs\n', blksz, t_psi);
 
end
tapp = toc(t00);
fprintf(' == Evalauting needlet approximation: J = %d, N = %d', J, N);
fprintf(', Time = %.2f secs\n', tapp);

sv = [sv_dir 'f' funtxt '_dN' '_QH' QH '_N' num2str(length(W_QH)) '_fis' num2str(kappa)...
            '_J' num2str(J) QNm '_L' num2str(deg_Nm) '.mat'];
save(sv,'dN')
end