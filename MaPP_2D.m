function [f, p, status] = MaPP_2D(T, n, K)

% [f, p, status] = MaPP_2D(T, n, K)
% 
% Matrix Pencil and Pairing (MaPP) carries out Vandermonde decomposition of a 2-level Toeplitz matrix T:
% T = sum_{k=1}^K p_k a(f_k)a^H(f_k)
% for any K <= N - max(n1,n2).
% 
% Inputs:
% T: 2-level Toeplitz matrix
% n=[n1, n2]: level dimensions
% K: # 2D sinusoids, or rank of T
% 
% Outputs:
% f: matrix of 2D frequencies, of dimension 2 * K
% p: power vector of frequencies
% status: indicator of success, 1: success, 0: failure
% 
% Reference:
% Z. Yang, L. Xie, and P. Stoica, "Multi-Dimensional Vandermonde Decomposition 
% and Its Use for Multi-Dimensional Super-Resolution", ISIT, 2015.
%
% Z. Yang, L. Xie and P. Stoica, "Vandermonde decomposition of multilevel Toeplitz 
% matrices with application to multidimensional super-resolution," 
% IEEE Transactions on Information Theory, 2016.
% 
% Written by Z. Yang, Jan 2015


status = true;

N = prod(n);

% eigendecomposition of T
[V0, Lambda] = eigs((T+T')/2, K);
V = V0 * sqrt(Lambda);

% f_1j
V1upp = V(1:N-n(2),:);
V1low = V(n(2)+1:end,:);
z1 = eig(V1upp'*V1low, V1upp'*V1upp);
f1 = sort(mod(imag(log(z1))/(2*pi), 1));

% f_2j
V2low = V;
V2upp = V;
V2low(1:n(2):end, :) = [];
V2upp(n(2):n(2):end, :) = [];
z2 = eig(V2upp'*V2low, V2upp'*V2upp);
f2 = sort(mod(imag(log(z2))/(2*pi), 1));

% pairing
idx = (1:K)'; % indices of f_2j for pairing with f_1j 
gval = zeros(K,1);
A1 = exp(1i*2*pi*kron((0:n(1)-1)',f1')) / sqrt(n(1));
A2 = exp(1i*2*pi*kron((0:n(2)-1)',f2')) / sqrt(n(2));
for j1 = 1:K
    for j2 = j1:K
        b = kron(A1(:,j1), A2(:,j2));
        gval(j2) = norm(b'*V0);
    end
    [m, kk] = max(gval(j1:K));
    kk = kk + j1 - 1;
    if m < .99
%         fprintf('Warning! low correlation.\n');
    end
    if kk ~= j1
        idx_temp = idx(j1);        
        idx(j1) = idx(kk);
        idx(kk) = idx_temp;
        
        a2_temp = A2(:,j1);
        A2(:,j1) = A2(:,kk);
        A2(:,kk) = a2_temp;
    end
end

mat_bbH = zeros(N^2,K);
for j = 1:K
    mat_bbH(:,j) = vec(kron(toeplitz(conj(A1(:,j))), toeplitz(conj(A2(:,j))))) / sqrt(N);
end
Tvec = T(:);
p = [real(mat_bbH); imag(mat_bbH)] \ [real(Tvec); imag(Tvec)];

res_rel = norm(Tvec - mat_bbH*p) / norm(Tvec);
if any(p<-1e-4) || res_rel>1e-6;
    status = false;
end

f = [f1 f2(idx)]';

end
