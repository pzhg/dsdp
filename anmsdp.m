function [f, varargout] = anmsdp(x, L, sigma, flag)
    % ANM SDP
    % By Zhe Zhang, 9/8/2016, zzhang18@gmu.edu
    %
    % Solve formulation: min ||x||_A. If there exists noise, the
    % formulation turns out to be: min ||x*||_A, s.t. ||x-x*||_2<=e
    %
    % The frequency components in x is given in f, which can be used to solve
    % DOA problems.
    %
    % x is a sinosoidal signal in the form of x=A*c, where A is array manifold
    % matrices and c is the amplitude of its frequency components.
    %
    % Need CVX installed, need MEMP_1D.m
    %
    % Input:
    %     x: Measurement
    %     L: Sparsity
    %     sigma: Noise Level. For noiseless case, let sigma=0.
    %     flag: DOA algorithm. 'music' or 'MEMP'.
    %
    % Output:
    %     f: Recovered (digital) frequencies / angles in range [0, 1]
    %     xe: Recovered (denoised) sinosoidal signal x. Only avaliable for 
    %        noise case.
    %     p: Recovered amplitide of sinosoidal components. Only avaliable
    %        for 'MEMP' method.
    
    N = length(x);

    %% Regularization Parameter
    if sigma > 0
        lambda = sigma * (1 + 1 / log(N)) * sqrt((N) * log(N) + (N) * log(4 * pi * log(N)));
    end

    %% SDP
    cvx_begin quiet
        
        variable v
        variable Tu(N, N) complex hermitian toeplitz
        variable xe(N, 1) complex

        if sigma > 0
            minimize lambda * ( v + 1 /  N * trace(Tu)) +  sum_square_abs(x - xe)
        else
            minimize v + 1 / N * trace(Tu)
        end

        subject to
            if sigma > 0
                [v, xe'; xe, Tu] == hermitian_semidefinite(N + 1)
            else
                [v, x'; x, Tu] == hermitian_semidefinite(N + 1)
            end

    cvx_end
    
    varargout = {};
    
    if sigma > 0
        varargout{end + 1} = xe;
    end
    
    %% DOA Estimation
    if strcmp(flag, 'music')
        f = mod(rootmusic(Tu, L) / 2 / pi, 1);
    elseif strcmp(flag, 'MEMP')
        [f, p, status] = MEMP_1D(Tu, N, L);
        varargout{end + 1} = p;
    end

end
