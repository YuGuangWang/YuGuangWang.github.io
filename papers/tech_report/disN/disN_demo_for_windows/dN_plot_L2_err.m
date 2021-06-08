function dN_plot_L2_err(J)
% plot maximum errors by contributions U_{jN} (discrete partial sum)
close all

ld_dir = 'disN_rbf\';
sv_dir_1 = ['dN' num2str(J)];
sv_dir = ['disN_fig\' sv_dir_1 '\'];
if ~exist(sv_dir,'dir')
    mkdir(sv_dir);
end
d = 2;
fis = 5;
% J = 6; %%%%%%
funtxt_ce = {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'};
QH = 'SD'; % discretisation quadrature %%%%%%
QNm = 'SD'; %%%%%% % evaluating quadrature
degJ = 3*2^(J-1) - 1;
deg_Nm = 301;
[W_QH,~] = pntset(QH,degJ);
[~,x_Nm] = pntset(QNm, deg_Nm);
fprintf('\n**************************************************************\n');
fprintf('Plot pointwise and L_2 errors for discrete needlet approximation\n');
fprintf('**************************************************************\n');

neord_j = 0:J;
fitn = 3; %%%%%%
nErr = 2.^neord_j;
nFit = 4:0.01:nErr(end);
err_L2 = zeros(length(funtxt_ce),length(neord_j));
err_L2_fit = zeros(length(funtxt_ce),length(nFit));
pstr_L2 = cell(1,length(neord_j));
lg_txt = cell(length(funtxt_ce),1);
for i_f = 1:length(funtxt_ce)
    funtxt = funtxt_ce{i_f};
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
    ld = [ld_dir 'f' funtxt '_dN' '_QH' QH '_N' num2str(length(W_QH))...
        '_fis' num2str(fis) '_J' num2str(J) QNm '_L' num2str(deg_Nm) '.mat'];
    % load discrete needlet approximation 'dN'
    load(ld, 'dN');
    % compute L2 errors
    f_Nm = fun(x_Nm',rbf_k)';
    switch funtxt
    case {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'}
        fun_fprt = 'normalised Wendland';
        k = rbf_k;
        delta = (3*k+3)*gamma(k+1/2)/(2*gamma(k+1));
    end
    N = size(x_Nm,2);
    dN_err = zeros(J+1,N);
    for j = 0:J
        dN_err(j+1,:) = abs(f_Nm - dN(j+1,:));
    end
    L2_err_dN = sqrt(sum(dN_err.^2,2)/N);
    % View for RBFs
    vw = [165 15];
    %%% Plot function    
    fprintf('\nTest function: %s, k = %d, delta = %.2f\n', fun_fprt, k, delta);
    s = k + 1.5;
    tstr_rbf = ['RBF f_' num2str(k) ' on S^2' 'with smoothness s=' num2str(s)];
    fprintf('***** Plot RBF f_%d *****\n', k);
    fig_rbf = pltfunc_rbf_loc_2(f_Nm,x_Nm);    
    view(vw)
    text(1.4,2,2.3,tstr_rbf);
    hold off
    sv_rbf = [sv_dir funtxt '_original'];
    % save fig_rbf to a jepg file
    print(fig_rbf,'-djpeg','-r300',sv_rbf);
    %%% Plot approxiamtion
    U = dN(end,:);
    fprintf('***** Plot discrete needlet approximation up to order J=%d *****\n', J);
    tstr_disN = ['Discrete needlet approxiamtion of order J=' num2str(J) ' for f_' num2str(k)];
    fig_disN = pltfunc_rbf_loc_2(U,x_Nm);
    text(1.7,2,2.3,tstr_disN);
    view(vw)
    hold off
    sv_disN = [sv_dir funtxt '_dN' num2str(J)];
    % save fig_rbf to a jepg file
    print(fig_disN,'-djpeg','-r300',sv_disN);
    % Calcualte and plot error
    E = f_Nm - U;
    abs_E = abs(E);
    Emax = max(max(abs_E));
    E2nrm = sqrt(sum(E.^2)/numel(E));
    fprintf(' === Approximation error for N = %d points: ', N);
    fprintf('Emax = %.2e, E2nrm = %.2e\n', Emax, E2nrm);
    fprintf('***** Plot discrete needlet approximation up to order J=%d *****\n', J);
    tstr_err = ['Errors of discrete needlet approxiamtion of order J=' num2str(J) ' for f_' num2str(k)];
    fig_err = pltfunc_rbf_loc_2(E,x_Nm);
    text(1.8,2,2.3,tstr_err);
    view(vw)
    sv_err = [sv_dir funtxt '_dN' num2str(J) '_err'];
    % save fig_rbf to a jepg file
    print(fig_err,'-djpeg','-r300',sv_err);
    % Curve fitting for L2 errors
    if J == 7
        switch funtxt
            case {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3'}
                [p, pstr] = fitpowJ(nErr(fitn:end),L2_err_dN(fitn:end));
            case {'nrWend_k4'}
                [p, pstr] = fitpowJ(nErr(fitn:end-1),L2_err_dN(fitn:end-1));
        end
    else
        [p, pstr] = fitpowJ(nErr(fitn:end),L2_err_dN(fitn:end));
    end
    Fit_err_L2 = p(1)*nFit.^p(2);
    err_L2(i_f,:) = L2_err_dN;
    err_L2_fit(i_f,:) = Fit_err_L2;
    pstr_L2{i_f} = pstr;
    % legend text
    lg_txt{i_f} = ['k = ' num2str(rbf_k) ', ' 's = ' num2str(rbf_k + 1.5)];
end
%%%%% Plot L2 error %%%%%
lw = 1.4;
mk = 8;
fig_L2err = figure;
set(0,'DefaultAxesColorOrder',...
    [0 0.498 0;0 0.498 0;
    1 0 0;1 0 0;...
    0 0.749 0.749;0 0.749 0.749;...
    0.498 0 0;0.498 0 0;
    0 0 0.498;0 0 0.498;
    0.749 0.749 0;0.749 0.749 0;
    ]);
loglog(nErr,err_L2(1,:),'v',...
    nFit,err_L2_fit(1,:),'-',...
    nErr,err_L2(2,:),'*',...
    nFit,err_L2_fit(2,:),'-',...
    nErr,err_L2(3,:),'p',...
    nFit,err_L2_fit(3,:),'-',...
    nErr,err_L2(4,:),'h',...
    nFit,err_L2_fit(4,:),'-',...
    nErr,err_L2(5,:),'o',...
    nFit,err_L2_fit(5,:),'-',...
    'LineWidth',lw,'MarkerSize',mk);
switch J
    case 4
        ylpow = -8;
    case 5
        ylpow = -11;
    case 6
        ylpow = -15;
end
ylminmax = [10^ylpow 10^0]; %%%%%%
ytick = 10.^(ylpow:2:0); %%%%%%
ax = gca;
set(gca,'XTick',[1 2 4 8 16 32 64 128]);
set(ax,'XTicklabel',{'0','1','2','3','4','5','6','7'});
xlabel('Order J');
ylabel('L_2 errors');
% title
ti = ['Discrete needlet approximation of order ' num2str(J) ' for RBFs on S^' num2str(d)];
title(ti);
xStar = 1.8; %%%%%%
xEnd = 2^(J+0.15);
xlim([xStar xEnd])
ylim(ylminmax) %%%%%%
set(gca,'YTick',ytick);
sv_funtxt = 'nrWk';
hold on
h_lg = legend(lg_txt{1},pstr_L2{1},...
    lg_txt{2},pstr_L2{2},...
    lg_txt{3},pstr_L2{3},...
    lg_txt{4},pstr_L2{4},...
    lg_txt{5},pstr_L2{5},...
    'Location','SouthWest');
set(h_lg,'FontSize',14);
grid on
set(gca,'YMinorGrid','off')
hold off
sv = [sv_dir 'f' sv_funtxt '_dN' '_errL2'...
    '_QH' QH '_N' num2str(length(W_QH)) '_QNm' QNm  '_L' num2str(deg_Nm)];
% save fig_L2err to a jepg file
fprintf('\n***** Fitting rates of convergence errors *****\n');
print(fig_L2err,'-djpeg','-r300',sv);