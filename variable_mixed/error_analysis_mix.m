function [max_error, convergence_rate] = error_analysis_mix(N_values)

    a = -1; b = 1;
    alpha_fun = @(x) 0.9 - 0.1*x.^2;
    %alpha_fun = @(x) 0.8 - 0.2*x.^2;
    
    v_fun = @(x) 0.5+0.2*x;
    mu = 1;
    u_exact = @(x) (1 - x.^2).^(alpha_fun(x) / 2);

    function d2F = d2u_exact(x)


        alpha = 0.9 - 0.1*x.^2;
        alpha_prime = -0.2*x;
        alpha_dprime = -0.2*ones(size(x));


       
        % alpha = 0.8 - 0.2*x.^2;
        % alpha_prime = -0.4*x;
        % alpha_dprime = -0.4*ones(size(x));
        

        tol = 1e-12;
        x_safe = x;
        x_safe(abs(x) >= 1 - tol) = sign(x(abs(x) >= 1 - tol)) * (1 - tol);
        G = (alpha_prime ./ 2) .* log(1 - x_safe.^2) - (alpha .* x_safe) ./ (1 - x_safe.^2);
        Gprime = 0.5 * alpha_dprime .* log(1 - x_safe.^2) ...
                 - (2 * alpha_prime .* x_safe + alpha) ./ (1 - x_safe.^2) ...
                 - (2 * alpha .* x_safe.^2) ./ ((1 - x_safe.^2).^2);
        F = (1 - x_safe.^2).^(alpha / 2);
        d2F = F .* (G.^2 + Gprime);
        d2F(abs(x) >= 1 - tol) = 0;
    end

    N_ref = 2^14;
    x_ref = graded_mesh(a, b, N_ref, alpha_fun);         
    x_ref_col = x_ref(:);
    x_ref_internal = x_ref(2:end-1)';                     

    A_ref = coef_matrix(N_ref, alpha_fun);                
    frac_full = A_ref * u_exact(x_ref_col);               
    frac_internal = frac_full(2:end-1);                   

    
    interp_frac = @(xq) interp1(x_ref_internal, frac_internal, xq(:), 'spline');

   
    max_error = zeros(1, length(N_values));
    convergence_rate = zeros(1, length(N_values)-1);

    for idx = 1:length(N_values)
        N = N_values(idx);
        x = graded_mesh(a, b, N, alpha_fun);
        x_internal = x(2:end-1)';                 

       
        A_full = coef_matrix_mixed(N, v_fun, mu, alpha_fun);
        F_full = zeros(N+1, 1);

        v = v_fun(x_internal);
        delta_u = d2u_exact(x_internal);
        frac_interp = interp_frac(x_internal);
        u_ex = u_exact(x_internal);
        f = -v .* delta_u + (1 - v) .* frac_interp + mu * u_ex;

        F_full(2:end-1) = f;

        u_h_full = A_full \ F_full;
        u_h_internal = u_h_full(2:end-1);
        error = max(abs(u_h_internal - u_ex));

        max_error(idx) = error;

        if idx > 1
            conv = log(max_error(idx-1) / max_error(idx)) / ...
                   log(N_values(idx) / N_values(idx-1));
            convergence_rate(idx-1) = conv;
        end
    end

    fprintf('max_error (混合方程):\n');
    fprintf(' %.2e', max_error);
    fprintf('\n');
    fprintf('convergence_rate:\n');
    fprintf(' %.2f', convergence_rate);
    fprintf('\n');
end
