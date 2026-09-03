
function x = graded_mesh(a, b, N, alpha_fun)

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
