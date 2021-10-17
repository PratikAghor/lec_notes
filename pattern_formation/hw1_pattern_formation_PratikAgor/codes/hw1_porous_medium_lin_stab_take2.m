%*************************************************************************%
% Linear stability of porous medium convection with salt
% Assuming z-variation to be sin(m*pi*z)
% Written by Pratik Aghor
%*************************************************************************%
close all; clear all; clc;
Ra_T=4*pi*pi+ 10;
Ra_S=0;
tau = 0.0; % diffusion coefficient for salt/ for temperature 

nk = 32;

kvec = linspace(-nk/2, nk/2, nk+1);
max_growth_rate_vs_k = zeros(nk+1, 2);
kx = 1; % for testing
%*************************************************************************%
%*************************************************************************%
%*************************************************************************%
m = 1;
for j = 1:nk+1
    kx = kvec(j);
    % order of perturbation variables: psi (streamfunction), theta (temperature), s (salt)
    ksq = kx*kx;
    msq_pisq = m*m*pi*pi;

    A11 = ksq + msq_pisq; A12 = -1i*kx*Ra_T; A13 = 1i*kx*Ra_S;
    A21 = -1i*kx; A22 = -(ksq + msq_pisq); A23 = 0;
    A31 = -1i*kx; A32 = 0; A33 = -tau*(ksq+ msq_pisq);
    
    %*************************************************************************%
    %build matrices
    A=[A11,A12,A13;A21,A22,A23;A31,A32,A33];
    B=[0, 0, 0; 0, 1, 0; 0, 0, 1];
    %*************************************************************************%
    eigvals=eig(A,B);
    % filter infinities
    idx = isfinite(eigvals);
    eigvals = sort(eigvals(idx), 'descend', 'ComparisonMethod', 'abs');
    max_growth_rate_vs_k(j, 1) = kx;
    max_growth_rate_vs_k(j, 2) = eigvals(1);
    %*************************************************************************%
end
%*************************************************************************%
figure(1)
plot(max_growth_rate_vs_k(:, 1), max_growth_rate_vs_k(:, 2), '-o', 'LineWidth', 2 )
hold on
plot(kvec, ((kvec.*kvec) ./(kvec.*kvec + msq_pisq))*Ra_T - (kvec.*kvec + msq_pisq), '-', 'LineWidth', 2)
legend({'2-components','1-component'},'Location','northEast')

xlabel('$k$','Interpreter','latex','FontSize',24); ylabel('$\sigma$','Interpreter','latex','FontSize',24);
grid on
title(['dispersion relation'])
%*************************************************************************%     