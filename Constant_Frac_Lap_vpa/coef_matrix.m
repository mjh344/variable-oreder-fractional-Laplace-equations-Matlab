function A = coef_matrix(N, alp)
% 使用VPA计算网格差和离散系数
% 最终将矩阵转成double，以便高效求解线性方程
%
% 要求：0 < alp < 1

    alp = vpa(alp);

    a = vpa(-1);
    b = vpa(1);

    % x是VPA符号数组
    x = graded_mesh(a, b, N, alp);

    % h也是VPA符号数组，不会因为端点舍入而出现0
    h = x(2:end)-x(1:end-1);

    c_alpha = alp * vpa(2)^(alp-1) ...
        * gamma((1+alp)/2) ...
        / (sqrt(sym(pi))*gamma(1-alp/2));

    d = 1/(alp*(1-alp));

    % 主对角线：先用VPA计算，再转成double
    aii_vpa = ...
        (x(2:end-1)-x(1:end-2)).^(1-alp)./h(1:end-1) ...
        +(x(3:end)-x(2:end-1)).^(1-alp)./h(2:end);

    aii = double(vpa(aii_vpa));

    % 上对角线
    aii1_vpa = ...
        -(h(3:end)+h(2:end-1)) ...
        .*(x(3:end-1)-x(2:end-2)).^(1-alp) ...
        ./h(2:end-1)./h(3:end) ...
        +(x(4:end)-x(2:end-2)).^(1-alp)./h(3:end);

    aii1 = double(vpa(aii1_vpa));

    % 下对角线
    aii_1_vpa = ...
        (x(3:end-1)-x(1:end-3)).^(1-alp)./h(1:end-2) ...
        -(h(1:end-2)+h(2:end-1)) ...
        .*(x(3:end-1)-x(2:end-2)).^(1-alp) ...
        ./h(1:end-2)./h(2:end-1);

    aii_1 = double(vpa(aii_1_vpa));

    temp_r = cell(N-1,1);
    temp_l = cell(N-1,1);
    ai = cell(N-1,1);
    ak = cell(N-1,1);

    % 右侧非三对角项
    for i = 1:N-3

        temp_vpa = ...
            (x(i+2:end-2)-x(i+1)).^(1-alp) ...
            ./h(i+2:end-1) ...
            -(h(i+3:end)+h(i+2:end-1)) ...
            .*(x(i+3:end-1)-x(i+1)).^(1-alp) ...
            ./(h(i+2:end-1).*h(i+3:end)) ...
            +(x(i+4:end)-x(i+1)).^(1-alp) ...
            ./h(i+3:end);

        % 每一行计算完成后再转成double
        temp_r{i} = double(vpa(temp_vpa));

        ai{i} = [zeros(1,i+1), temp_r{i}];
    end

    % 左侧非三对角项
    for k = 3:N-1

        temp_vpa = ...
            (x(k+1)-x(1:k-2)).^(1-alp)./h(1:k-2) ...
            -(h(2:k-1)+h(1:k-2)) ...
            .*(x(k+1)-x(2:k-1)).^(1-alp) ...
            ./(h(1:k-2).*h(2:k-1)) ...
            +(x(k+1)-x(3:k)).^(1-alp)./h(2:k-1);

        temp_l{k} = double(vpa(temp_vpa));

        ak{k} = [temp_l{k}, zeros(1,N-k+1)];
    end

    % 组装双精度矩阵
    A1 = diag(aii) ...
        +diag(aii_1,-1) ...
        +diag(aii1,1) ...
        +vertcat(ai{1:N-3}, zeros(2,N-1)) ...
        +vertcat(zeros(2,N-1), ak{3:N-1});

    % 归一化系数也先用VPA计算
    factor = double(vpa(c_alpha*d));

    A1 = factor*A1;

    if any(~isfinite(A1(:)))
        error('矩阵中仍然存在Inf或NaN，请增加digits设置。');
    end

    % 齐次Dirichlet边界条件
    A = [zeros(N-1,1), A1, zeros(N-1,1)];

    A = [1, zeros(1,N);
         A;
         zeros(1,N), 1];
end