function [f, H]=d2dsdp(Y, XR, XT, L, N, M, Tr, Tt, sigma, flag)

% Decoupled 2D SDP
% By Zhe Zhang, 9/8/2016, zzhang18@gmu.edu
% Solve formulation: Y=H*X, need CVX installed, need MEMP_1D.m
% Input:
%     Y: Measurement, M-by-T
%     X: Training set, N-by-T. For DOA, let X=eye(N)
%     L: Sparsity
%     N, M: Size of H
%     T: Size of Training. For DOA, let T=N
%     sigma: Noise Level. For noiseless case, let sigma=0
%     flag: DOA algorithm. 'music' or 'MEMP'.
% Output:
%     f: Frequencies
%     H: Channel

%% Regularization Parameter
if sigma>0
    lambda=sigma*(1+1/log(M+N))*sqrt((M+N)*log(M+N)+(M+N)*log(4*pi*log(M+N)));
end

%% Decoupled 2D SDP
cvx_begin quiet
    variable H(N, M) complex
    variable Tu_r(N, N) complex hermitian toeplitz
    variable Tu_t(M, M) complex hermitian toeplitz
    if sigma>0
        minimize lambda/(2*sqrt(M*N))*(trace(Tu_r)+trace(Tu_t))+1/2*quad_form(Y(:)-reshape(XR*H*XT, [Tr*Tt, 1]), eye(Tr*Tt))
    else
        minimize 1/(2*sqrt(M*N))*(trace(Tu_r)+trace(Tu_t))
    end
    subject to
        if sigma>0
            [Tu_t, H'; H, Tu_r]==hermitian_semidefinite(M+N)
        else
            [Tu_t, H'; H, Tu_r]==hermitian_semidefinite(M+N)
            Y==H*X
        end
cvx_end

%% DOA Estimation
if strcmp(flag, 'music')
    fr_e=mod(rootmusic(Tu_r, L)/2/pi, 1);
    ft_e=mod(rootmusic(Tu_t, L)/2/pi, 1);
elseif strcmp(flag, 'MEMP')
    [ft_e, p, status]=MEMP_1D(Tu_t, M, L);
    [fr_e, p, status]=MEMP_1D(Tu_r, N, L);
end
fr_e=fr_e';
ft_e=ft_e';
f=[fr_e; ft_e];

%% Pairing
v_N=[0:(N-1)]';
v_M=[0:(M-1)]';
A_r=[];
for index=1:L
    A_r=[A_r, exp(1i*2*pi*fr_e(index)*v_N)];
end
Dr=diag(diag(abs(pinv(A_r)*Tu_r*pinv(A_r)')));
V=inv(Dr)*pinv(A_r)*H;
fr_sort=zeros(1, L);
for index=1:L
    c=abs(V*exp(1i*2*pi*ft_e(index)*v_M));
    [m, ii]=max(c);
    fr_sort(index)=fr_e(ii);
end
f=[fr_sort; ft_e];