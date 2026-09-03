clear; clc; close all;

% ---------- 参数设置 ----------
a = -1; b = 1;
N = 512;                          % 待测试网格点数
N_ref = 2^14;                     % 精细参考网格

% ---------- 定义变阶函数 ----------
alpha_funs = {
    %@(x) 0.9 - 0.1*abs(x),               % 线性
    @(x) 0.8 - 0.1*abs(x).^2,            % 二次
    @(x) 0.9 - 0.1*sin(abs(pi*x))        % 正弦
};
legend_str = {
              '$\alpha(x)=0.8-0.1x^2$', ...
              '$\alpha(x)=0.9-0.1\sin(\pi|x|)$'};
colors = {'r', 'b'};
line_styles = {'-', '--'};

% 创建图形
figure('Position', [100, 100, 600, 400]);
hold on;

% 对每个 alpha 计算并绘制误差
for k = 1:length(alpha_funs)
    alpha_fun = alpha_funs{k};
    
    % ---- 目标网格（使用依赖 alpha 的 r） ----
    x = graded_mesh(a, b, N, alpha_fun);
    x = x(:);
    A = coef_matrix(N, alpha_fun);
    
    % ---- 参考网格（同样使用依赖 alpha 的 r） ----
    x_ref = graded_mesh(a, b, N_ref, alpha_fun);
    x_ref = x_ref(:);
    A_ref = coef_matrix(N_ref, alpha_fun);
    u_exact_ref = (1 - x_ref.^2).^(alpha_fun(x_ref) / 2);
    f_ref = A_ref * u_exact_ref;
    
    % ---- 插值到目标网格 ----
    f_interp = @(xq) interp1(x_ref, f_ref, xq(:), 'spline');
    f = f_interp(x);
    
    % ---- 求解并计算误差 ----
    u_h = A \ f;
    u_exact = (1 - x.^2).^(alpha_fun(x) / 2);
    abs_error = abs(u_h - u_exact);
    
    % 提取内部点（边界点误差为零）
    x_internal = x(2:end-1);
    error_internal = abs_error(2:end-1);
    
    % ---- 绘制误差分布 ----
    semilogy(x_internal, error_internal, ...
             'Color', colors{k}, 'LineStyle', line_styles{k}, ...
             'LineWidth', 1.8, 'DisplayName', legend_str{k});
end

% ---------- 图形修饰 ----------
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$|u_h - u|$', 'Interpreter', 'latex', 'FontSize', 14);
title(sprintf('误差分布对比 (N=%d, r=2(2-\\alpha_{max}/\\alpha_{min}))', N), 'FontSize', 14);
legend('Interpreter', 'latex', 'Location', 'best');
grid on;
hold off;