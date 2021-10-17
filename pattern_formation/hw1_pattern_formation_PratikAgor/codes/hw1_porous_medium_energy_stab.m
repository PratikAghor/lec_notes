%*************************************************************************%
% Energy stability of porous medium convection with salt
% Assuming z-variation to be sin(m*pi*z)
% Written by Pratik Aghor
%*************************************************************************%
close all; clear all; clc;
Ra_T=4*pi*pi;
Ra_S=0;
tau = 0.0; % diffusion coefficient for salt/ for temperature 

nkmax = 16;
nk = 100;
kvec = linspace(-nkmax/2, nkmax/2, nk);

Lambda_vs_k = zeros(nk, 2);
%*************************************************************************%
m = 1;
for j = 1:nk
    kx = kvec(j);

    ksq = kx*kx;
    msq_pisq = m*m*pi*pi;
    lap_op = -(ksq + msq_pisq); % negative of discrete laplacian operator

    % order of variables: theta_E, c_E, v_E, w_E: check latex hw for
    % details

    A11 = -2*lap_op; A12 = 0; A13 = -ksq*Ra_T; A14 = -1;
    A21 = 0; A22 = -2*tau*lap_op; A23 = ksq*Ra_S; A24 = -1;
    A31 = 1; A32 = 1; A33 = lap_op; A34 = 0;
    A41 = ksq*Ra_T; A42 = -ksq*Ra_S; A43 = 0; A44 = lap_op;
    
    %*************************************************************************%
    %build matrices
    A=[A11, A12, A13, A14; A21, A22, A23, A24; A31, A32, A33, A34; A41, A42, A43, A44];
    B=[1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0];
    %*************************************************************************%
    eigvals=eig(A,B);
    idx = isfinite(eigvals);
    eigvals = sort(eigvals(idx), 'descend', 'ComparisonMethod', 'abs');
    Lambda_vs_k(j, 1) = kx;
    Lambda_vs_k(j, 2) = eigvals(1);
end

%*************************************************************************%
figure(1)
plot(Lambda_vs_k(:, 1), Lambda_vs_k(:, 2), '-o', 'LineWidth', 2 )
hold on
% legend({'2-components','1-component'},'Location','northEast')

xlabel('$k$','Interpreter','latex','FontSize',24); ylabel('$\Lambda$','Interpreter','latex','FontSize',24);
grid on
title(['energy stability curve'])
%*************************************************************************%     