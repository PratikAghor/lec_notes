%*************************************************************************%
% Linear stability of porous medium convection with salt
% This one does not assume z-variation to be sim(m*pi*z)
% Written by Pratik Aghor
%*************************************************************************%
close all; clear all; clc;
Ra_T=4*pi*pi;
Ra_S=0;
tau = 0; % diffusion coefficient for salt/ for temperature 

Nz=32;
nz = Nz + 1; % size of everything full
n = Nz - 1; % size of everything after cropping first and last row-cols

nkmax = 16;
nk = 100;
[D,y]=cheb(Nz);

% co-ordinate transform: y = a*z + b
a = 2; b = -1; 
z = (y-b)/a;

% chain rule d/dz = (d/dy)*(dy/dz) = a*D
D = a*D;
D2 = D*D;

% D = D(2:Nz, 2:Nz);
% D2 = D2(2:Nz, 2:Nz);

kvec = linspace(-nkmax/2, nkmax/2, nk);
max_growth_rate_vs_k = zeros(nk, 2);
% kx = 1; % for testing
%*************************************************************************%
%*************************************************************************%
Z=zeros(nz, nz); 

A11=zeros(nz, nz); A12=zeros(nz, nz); A13=zeros(nz, nz);
A21=zeros(nz, nz); A22=zeros(nz, nz); A23=zeros(nz, nz);
A31=zeros(nz, nz); A32=zeros(nz, nz); A33=zeros(nz, nz);
% B11=zeros(nz, nz); B12=zeros(nz, nz); B13=zeros(nz, nz);
% B21=zeros(nz, nz); B22=zeros(nz, nz); B23=zeros(nz, nz);
% B31=zeros(nz, nz); B32=zeros(nz, nz); B33=zeros(nz, nz);
%*************************************************************************%
for i = 1:nk
    kx = kvec(i);
    % order of perturbation variables: psi (streamfunction), theta (temperature), s (salt)
    ksq = kx*kx;
    A11 =  ksq*eye(nz) - D2; A12 = -1j*kx*Ra_T*eye(nz); A13 = 1j*kx*Ra_S*eye(nz);
    A21 = -1j*kx*eye(nz); A22 = -(ksq*eye(nz) - D2); 
    A31 = -1j*kx*eye(nz); A33 = -tau*(ksq*eye(nz)-D2);
    
    B22 = eye(nz);
    B33 = eye(nz); 
    %*************************************************************************%
    %build matrices
    A=[A11,A12,A13;A21,A22,A23;A31,A32,A33];
    B=[Z,Z,Z;Z,B22,Z;Z,Z,B33];

    %BCs
    A(1, :) = 0; A(1, 1) = 1; B(1, :) = 0;
    A(nz, :) = 0; A(nz, nz) = 1; B(nz, :) = 0;

    A(nz+1, :) = 0; A(nz+1, nz+1) = 1; B(nz+1, :) = 0;
    A(2*nz, :) = 0; A(2*nz, 2*nz) = 1; B(2*nz, :) = 0;

    A(2*nz+1, :) = 0; A(2*nz+1, 2*nz+1) = 1; B(2*nz+1, :) = 0;
    A(3*nz, :) = 0; A(3*nz, 3*nz) = 1; B(3*nz, :) = 0;
    %*************************************************************************%
    eigvals=eig(A,B);
    % filter infinities
    idx = isfinite(eigvals);
    eigvals = real(eigvals(idx)) ; % sort(eigvals(idx), 'descend');
    eigvals = sort(eigvals, 'descend');
    max_growth_rate_vs_k(i, 1) = kx;
    if(max(eigvals) ==0 )
        max_growth_rate_vs_k(i, 2) = eigvals(nz-1);
    else
        max_growth_rate_vs_k(i, 2) = max(eigvals);
    end
    %*************************************************************************%
end
%*************************************************************************%
figure(1)
plot(max_growth_rate_vs_k(:, 1), max_growth_rate_vs_k(:, 2), '-o')
xlabel('$k$','Interpreter','latex','FontSize',24); ylabel('$\sigma$','Interpreter','latex','FontSize',24);
grid on
title(['dispersion relation'])
%*************************************************************************%     