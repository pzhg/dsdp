clc;
clear all;
close all;
warning off;

%% Set Parameters
N = 16; % Size N
M = 16; % Size M
L = 4; % Sparsity (rank)
SNR = 20;

%% Frequencies
f_r = [0.35 0.1 0.67 0.92];
f_t = [0.55 0.34 0.87 0.06];
c = [12, 8, 10, 11];

%% Generate Array Manifolds
A_r = [];
A_t = [];
v_M = [0:(M - 1)]';
v_N = [0:(N - 1)]';

for ii = 1:L
    A_r = [A_r, exp(1i * 2 * pi * f_r(ii) * v_N)];
    A_t = [A_t, exp(1i * 2 * pi * f_t(ii) * v_M)];
end

H = A_r * diag(c) * A_t';

HW = awgn(H, SNR, 'measured');
W = HW - H;
sigma = sqrt(sum(abs(W(:)).^2) / length(W(:)));

%% Decoupled ANM
tic; [f, H_e] = d2dsdp(HW, eye(N), eye(M), L, N, M, N, M, sigma, 'music'); toc;
ef_r = f(1, :)
ef_t = f(2, :)

%% Vectorized ANM
tic; [f, H_e] = vsdp(HW, eye(N), eye(M), L, N, M, N, M, sigma); toc;
ef_r = f(1, :)
ef_t = f(2, :)
