function A = coef_matrix(N, alp)

    alp = vpa(alp);

    a = vpa(-1);
    b = vpa(1);

    x = graded_mesh(a, b, N, alp);

    
    h = x(2:end)-x(1:end-1);

    c_alpha = alp * vpa(2)^(alp-1) ...
        * gamma((1+alp)/2) ...
        / (sqrt(sym(pi))*gamma(1-alp/2));

    d = 1/(alp*(1-alp));

    
    aii_vpa = ...
        (x(2:end-1)-x(1:end-2)).^(1-alp)./h(1:end-1) ...
        +(x(3:end)-x(2:end-1)).^(1-alp)./h(2:end);

    aii = double(vpa(aii_vpa));

    
    aii1_vpa = ...
        -(h(3:end)+h(2:end-1)) ...
        .*(x(3:end-1)-x(2:end-2)).^(1-alp) ...
        ./h(2:end-1)./h(3:end) ...
        +(x(4:end)-x(2:end-2)).^(1-alp)./h(3:end);

    aii1 = double(vpa(aii1_vpa));

    
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

    
    for i = 1:N-3

        temp_vpa = ...
            (x(i+2:end-2)-x(i+1)).^(1-alp) ...
            ./h(i+2:end-1) ...
            -(h(i+3:end)+h(i+2:end-1)) ...
            .*(x(i+3:end-1)-x(i+1)).^(1-alp) ...
            ./(h(i+2:end-1).*h(i+3:end)) ...
            +(x(i+4:end)-x(i+1)).^(1-alp) ...
            ./h(i+3:end);

        
        temp_r{i} = double(vpa(temp_vpa));

        ai{i} = [zeros(1,i+1), temp_r{i}];
    end

    
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

    A1 = diag(aii) ...
        +diag(aii_1,-1) ...
        +diag(aii1,1) ...
        +vertcat(ai{1:N-3}, zeros(2,N-1)) ...
        +vertcat(zeros(2,N-1), ak{3:N-1});

    factor = double(vpa(c_alpha*d));

    A1 = factor*A1;

    A = [zeros(N-1,1), A1, zeros(N-1,1)];

    A = [1, zeros(1,N);
         A;
         zeros(1,N), 1];
end
