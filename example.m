clc;
clear all;
close all;
warning off;

N=10;
M=10;
L=4;
SNR=20;

f_r=[0.35 0.1 0.67 0.92];
f_t=[0.55 0.34 0.87 0.06];
c=[12, 8, 10, 11];

A_r=[];
A_t=[];
v_M=[0:(M-1)]';
v_N=[0:(N-1)]';
for ii=1:L
    A_r=[A_r, exp(1i*2*pi*f_r(ii)*v_N)];
    A_t=[A_t, exp(1i*2*pi*f_t(ii)*v_M)];
end
H=A_r*diag(c)*A_t';

HW=awgn(H, SNR, 'measured');
W=HW-H;
sigma=sqrt(sum(abs(W(:)).^2)/length(W(:)));
[f, H_e]=d2dsdp(HW, eye(N), eye(M), L, N, M, N, M, sigma, 'music');

ef_r=f(1, :)
ef_t=f(2, :)
