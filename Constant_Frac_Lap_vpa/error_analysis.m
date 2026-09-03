function [max_error, convergence_rate] = ...
    error_analysis(N_values, alp_values)

    max_error = zeros(length(alp_values), length(N_values));

    convergence_rate = zeros( ...
        length(alp_values), length(N_values)-1);

    a = vpa(-1);
    b = vpa(1);

    for alp_idx = 1:length(alp_values)

        alp = vpa(alp_values(alp_idx));

        exact_constant = sqrt(sym(pi)) ...
            /(vpa(2)^alp ...
            *gamma(1+alp/2) ...
            *gamma((1+alp)/2));

        for N_idx = 1:length(N_values)

            N = N_values(N_idx);

            fprintf('正在计算：alp = %.2f, N = %d\n', ...
                double(alp), N);

            
            A = coef_matrix(N, alp);

            % 右端项
            f = zeros(N+1,1);
            f(2:N) = 1;

            u = A\f;

            x_vpa = graded_mesh(a, b, N, alp);

            one_minus_x2 = 1-x_vpa.^2;

            u_exact_vpa = exact_constant ...
                .*one_minus_x2.^(alp/2);

            u_exact = double(vpa(u_exact_vpa(:)));

            error_value = max(abs(u-u_exact));

            max_error(alp_idx,N_idx) = error_value;

            if N_idx > 1
                convergence_rate(alp_idx,N_idx-1) = ...
                    log(max_error(alp_idx,N_idx-1) ...
                    /max_error(alp_idx,N_idx)) ...
                    /log(N_values(N_idx) ...
                    /N_values(N_idx-1));
            end
        end
    end

    disp('max_error:');

    for alp_idx = 1:length(alp_values)
        fprintf('alp = %.1f:', double(alp_values(alp_idx)));
        fprintf(' %.6e', max_error(alp_idx,:));
        fprintf('\n');
    end

    disp('convergence_rate:');

    for alp_idx = 1:length(alp_values)
        fprintf('alp = %.1f:', double(alp_values(alp_idx)));
        fprintf(' %.4f', convergence_rate(alp_idx,:));
        fprintf('\n');
    end
end
