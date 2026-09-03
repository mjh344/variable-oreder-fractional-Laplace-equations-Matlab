
% 构造混合扩散方程的离散矩阵（内部节点）
%   -v(x) Δu + (1-v(x))(-Δ)^{α/2}u + μu = f,   u(±1)=0
% 输入：
%   N         : 网格分段数（节点数为 N+1，内部节点数为 N-1）
%   alpha_fun : 阶函数句柄，如 @(x) 0.9-0.1*x.^2
%   v_fun     : 速度函数句柄，默认 @(x) x
%   mu        : 反应系数，默认 1
% 输出：
%   A         : (N-1)×(N-1) 稀疏矩阵
function A = coef_matrix_mixed(N, v_fun, mu, alpha_fun)

    if nargin < 4
        alpha_fun = @(x) 0.9 - 0.1*x.^2;
    end

    a = -1; b = 1;
    x = graded_mesh(a, b, N, alpha_fun);   % 行向量，长度 N+1
    h = x(2:end) - x(1:end-1);             % 步长
    M = N - 1;                             % 内部节点数
    x_internal = x(2:end-1)';              % 内部节点列向量

    % ---------- 构造二阶导数矩阵 L_diff （内部节点） ----------
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

    % ---------- 提取分数阶内部矩阵 ----------
    A_full_frac = coef_matrix(N, alpha_fun);        % (N+1)×(N+1)
    A_internal = A_full_frac(2:end-1, 2:end-1);     % M×M

    % ---------- 组合内部矩阵 ----------
    v = v_fun(x_internal);
    A_internal_mixed = diag(-v) * L_diff + diag(1 - v) * A_internal + mu * speye(M);

    % ---------- 扩展为完整矩阵（边界条件） ----------
    A = sparse(N+1, N+1);
    A(2:end-1, 2:end-1) = A_internal_mixed;
    A(1, 1) = 1;
    A(end, end) = 1;
end