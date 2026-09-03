function [max_error, convergence_rate] = error_analysis_ref(N_values,alpha_fun)
% 使用参考解（N=4096）代替精确解，计算不同 N 下的最大误差和收敛阶
    a = -1;
    b = 1;

    % ---------- 1. 计算参考解 (N_ref = 2^14) ----------
    N_ref = 2^14;
    x_ref = graded_mesh(a, b, N_ref,alpha_fun);                     % 参考解网格
    
    A_ref = coef_matrix(N_ref);                    % 参考解系数矩阵
    f_ref = zeros(N_ref+1, 1);
    f_ref(2:end-1) = sin(x_ref(2:end-1));
    %f_ref(2:end-1) = x_ref(2:end-1);               % 右端项为 x
    %f_ref(2:end-1) = x_ref(2:end-1).^2;            % 右端项为 x^2
    


    u_ref = A_ref \ f_ref;                         % 参考数值解

    %构造插值函数（样条插值，光滑解下精度较高）
    u_ref_interp = @(xq) interp1(x_ref, u_ref, xq, 'spline');

    % 修改为分段线性插值
    %u_ref_interp = @(xq) interp1(x_ref, u_ref, xq, 'linear');

    % ---------- 2. 对每个 N 计算误差和收敛阶 ----------
    max_error = zeros(1, length(N_values));
    convergence_rate = zeros(1, length(N_values)-1);

    for N_idx = 1:length(N_values)
        N = N_values(N_idx);
                               
        x = graded_mesh(a, b, N, alpha_fun);
        A = coef_matrix(N);                        % 当前系数矩阵
        f = zeros(N+1, 1);
        f(2:end-1) = sin(x(2:end-1));
        %f(2:end-1) = x(2:end-1);                   % 右端项为 x
        %f(2:end-1) = x(2:end-1).^2;                % 右端项为 x^2
        
        
        u = A \ f;                                 % 当前数值解

        % 参考解在当前网格点上的值
        u_exact_at_x = u_ref_interp(x);
        u_exact_at_x = u_exact_at_x(:); 

        % 最大模误差
        error = max(abs(u - u_exact_at_x));
        max_error(N_idx) = error;

        % 计算收敛阶（对数比值）
        if N_idx > 1
            convergence_rate(N_idx-1) = log(max_error(N_idx-1) / max_error(N_idx)) / ...
                                       log(N_values(N_idx) / N_values(N_idx-1));
        end
    end

    % 输出结果
    disp('max_error:');
    fprintf(' %.2e', max_error);
    fprintf('\n');
    disp('convergence_rate:');
    fprintf(' %.2f', convergence_rate);
    fprintf('\n');
end