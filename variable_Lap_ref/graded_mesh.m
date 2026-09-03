
function x = graded_mesh(a, b, N, alpha_fun)
% 生成非均匀网格（基于变阶函数 alpha）
% 输入：a,b 区间端点，N 网格数（内部节点数），alpha_fun 变阶函数句柄
% 输出：网格点向量（长度 N+1）

    % 计算 alpha 在 [a,b] 上的最大值和最小值
    % alpha_min = fminbnd(alpha_fun, a, b);
    % alpha_max = -fminbnd(@(x) -alpha_fun(x), a, b);

    x_sample = linspace(a, b, 100001);
    alpha_values = alpha_fun(x_sample);

    alpha_min = min(alpha_values);
    alpha_max = max(alpha_values);

    %r = 2 * (2 - alpha_max) / alpha_min;
    r =1;
    
    
    j = 0:N;
    x = zeros(size(j));
    left_indices = j <= N/2;
    x(left_indices) = a + (b - a)/2 * (j(left_indices)/(N/2)).^r;
    x(~left_indices) = b - (b - a)/2 * (2 - j(~left_indices)/(N/2)).^r;
end