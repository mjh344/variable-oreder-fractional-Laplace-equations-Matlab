
function A = coef_matrix_mixed(N, v_fun, mu, alpha_fun)

    if nargin < 4
        alpha_fun = @(x) 0.9 - 0.1*x.^2;
    end

    a = -1; b = 1;
    x = graded_mesh(a, b, N, alpha_fun);   
    h = x(2:end) - x(1:end-1);            
    M = N - 1;                             
    x_internal = x(2:end-1)';              

    L_diff = sparse(M, M);
    for i = 1:M
        hL = h(i);
        hR = h(i+1);
        denom = hL + hR;
        if i > 1
            L_diff(i, i-1) = 2 / (denom * hL);
        end
        L_diff(i, i) = -2 / denom * (1/hL + 1/hR);
        if i < M
            L_diff(i, i+1) = 2 / (denom * hR);
        end
    end

   
    A_full_frac = coef_matrix(N, alpha_fun);        
    A_internal = A_full_frac(2:end-1, 2:end-1);     

    
    v = v_fun(x_internal);
    A_internal_mixed = diag(-v) * L_diff + diag(1 - v) * A_internal + mu * speye(M);

    A = sparse(N+1, N+1);
    A(2:end-1, 2:end-1) = A_internal_mixed;
    A(1, 1) = 1;
    A(end, end) = 1;
end
