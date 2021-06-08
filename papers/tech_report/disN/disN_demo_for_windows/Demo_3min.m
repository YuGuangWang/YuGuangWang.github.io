% Demo for discrete needlet approximations on S^d for RBFs
clear,clc
close all

% compute discrete needlet approximation up to level 5
J = 5; % order J
NeedletSBF(J);

% plot figures
dN_plot_L2_err(J);

fprintf('\n***** Algorithm ends! *****\n');