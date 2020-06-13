function [f, H]=vsdp(Y, XR, XT, L, N, M, Tr, Tt, sigma)

% Vectorized 2D SDP
% By Zhe Zhang, 9/8/2016, zzhang18@gmu.edu
% Solve formulation: Y=XR*H*XT', need CVX installed, need MaPP_2D.m
% Input:
%     Y: Measurement, Tr-by-Tt
%     XR: Receiver Modifier, Tr-by-N. For DOA, let XR=eye(N)
%	  XT: Transmitter Modifier, Tt-by-M. For DOA, let XT=eye(M)
%     L: Sparsity
%     N, M: Size of H
%     Tr: Size of XR. For DOA, let Tr=N
%     Tt: Size of XT. For DOA, let Tt=M
%     sigma: Noise Level. For noiseless case, let sigma=0
% Output:
%     f: Frequencies
%     H: Channel

%% Regularization Parameter
if sigma>0
    lambda=sigma*(1+1/log(M*N))*sqrt((M*N)*log(M*N)+(M*N)*log(4*pi*log(M*N)));
end

%% Vectorized 2D SDP
y=Y(:);
XX=kron(XT.', XR);
cvx_begin quiet
    variable v
    variable h(N*M, 1) complex
    variable Tu(N*M, N*M) complex hermitian
    variable u(2*N-1, M) complex
    if sigma>0
        minimize lambda/2*(trace(Tu)+v)+1/2*quad_form(y-XX*h, eye(Tr*Tt))
    else
        minimize 1/2*(trace(Tu)+v)
    end
    subject to
        if sigma>0
            [v, h'; h, Tu]==hermitian_semidefinite(N*M+1)
            for nt=1:M
                for nr=nt:M
                    Tu((N*(nt-1)+1):(N*nt), (N*(nr-1)+1):(N*nr))==toeplitz(u(N:end, nr-nt+1), u(N:-1:1, nr-nt+1));
                end
            end
        else
            [v, h'; h, Tu]==hermitian_semidefinite(N*M+1)
            for nt=1:M
                for nr=nt:M
                    Tu((N*(nt-1)+1):(N*nt), (N*(nr-1)+1):(N*nr))==toeplitz(u(N:end, nr-nt+1), u(N:-1:1, nr-nt+1));
                end
            end
            y==XX*h
        end
cvx_end

%% DOA Estimation and Pairing
[ff, p, status]=MaPP_2D(Tu, [N, M], L);
fr_e=ff(2, 1:L);
ft_e=1-ff(1, 1:L);
f=[fr_e; ft_e];
H=reshape(h, [N, M]);
