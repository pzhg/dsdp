function [f, p, status] = MEMP_1D(T, N, K)

    % Matrix Pencil from Zai Yang's Paper
    % By Zhe Zhang
    % T: Toeplitz Matrix; N: Matrix Size; K: Matrix Rank

    status = true;

    [V0, Lambda] = eigs((T + T') / 2, K);
    V = V0 * sqrt(Lambda);

    V1upp = V(1:N - 1, :);
    V1low = V(2:end, :);
    z1 = eig(V1upp' * V1low, V1upp' * V1upp);
    f = sort(mod(imag(log(z1)) / (2 * pi), 1));

    A1 = exp(1i * 2 * pi * kron((0:N - 1)',f')) / sqrt(N);
    mat_bbH = zeros(N^2, K);

    for j = 1:K
        mat_bbH(:, j) = vec(toeplitz(conj(A1(:, j)))) / sqrt(N);
    end

    Tvec = T(:);
    p = [real(mat_bbH); imag(mat_bbH)] \ [real(Tvec); imag(Tvec)];

    res_rel = norm(Tvec - mat_bbH * p) / norm(Tvec);

    if any(p <- 1e-4) || res_rel > 1e-6;
        status = false;
    end

end