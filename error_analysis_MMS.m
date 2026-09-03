function [max_error, convergence_rate] = error_analysis_MMS(N_values,alpha_fun)
% 制造解方法，精确解 u(x) = (1 - x^2)^(alpha(x)/2)
    a = -1; b = 1;

    % ---------- 变阶函数与精确解 ----------
    
    %alpha_fun = @(x) 0.9 - 0.1*sin(abs(pi*x));
    %alpha_fun = @(x) 0.8 - 0.1*abs(x).^2;
    u_exact = @(x) (1 - x.^2).^(alpha_fun(x) / 2);

    % ---------- 在参考网格上计算 f = L u ----------
    N_ref = 2^14;
    x_ref = graded_mesh(a, b, N_ref,alpha_fun);
    x_ref = x_ref(:);                     
    A_ref = coef_matrix(N_ref);
    u_ref_exact = u_exact(x_ref);        
    f_ref = A_ref * u_ref_exact;          

    % 插值函数（查询点强制列向量，返回列向量）
    f_interp = @(xq) interp1(x_ref, f_ref, xq(:), 'spline');

    % ---------- 对每个待测 N 求解 ----------
    max_error = zeros(1, length(N_values));
    convergence_rate = zeros(1, length(N_values)-1);

    for idx = 1:length(N_values)
        N = N_values(idx);
        x = graded_mesh(a, b, N,alpha_fun);
        x = x(:);                         
        A = coef_matrix(N);
        f = f_interp(x);                 
        u_h = A \ f;                      
        u_ex = u_exact(x);                

        error = max(abs(u_h - u_ex));     
        max_error(idx) = error;

        if idx > 1
            conv = log(max_error(idx-1) / max_error(idx)) / ...
                   log(N_values(idx) / N_values(idx-1));
            convergence_rate(idx-1) = conv;
        end
    end

    % 输出
    disp('max_error:');
    fprintf(' %.2e', max_error);
    fprintf('\n');
    disp('convergence_rate:');
    fprintf(' %.2f', convergence_rate);
    fprintf('\n');
end