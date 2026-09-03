function [max_error, convergence_rate] = error_analysis_ref(N_values,alpha_fun)

    a = -1;
    b = 1;

   
    N_ref = 2^14;
    x_ref = graded_mesh(a, b, N_ref,alpha_fun);                    
    
    A_ref = coef_matrix(N_ref);                    
    f_ref = zeros(N_ref+1, 1);
    f_ref(2:end-1) = sin(x_ref(2:end-1));
        
    u_ref = A_ref \ f_ref;                         

    
    u_ref_interp = @(xq) interp1(x_ref, u_ref, xq, 'spline');

    %u_ref_interp = @(xq) interp1(x_ref, u_ref, xq, 'linear');

   
    max_error = zeros(1, length(N_values));
    convergence_rate = zeros(1, length(N_values)-1);

    for N_idx = 1:length(N_values)
        N = N_values(N_idx);
                               
        x = graded_mesh(a, b, N, alpha_fun);
        A = coef_matrix(N);                      
        f = zeros(N+1, 1);
        f(2:end-1) = sin(x(2:end-1));
       
        u = A \ f;                                

        
        u_exact_at_x = u_ref_interp(x);
        u_exact_at_x = u_exact_at_x(:); 

        
        error = max(abs(u - u_exact_at_x));
        max_error(N_idx) = error;

        if N_idx > 1
            convergence_rate(N_idx-1) = log(max_error(N_idx-1) / max_error(N_idx)) / ...
                                       log(N_values(N_idx) / N_values(N_idx-1));
        end
    end

    disp('max_error:');
    fprintf(' %.2e', max_error);
    fprintf('\n');
    disp('convergence_rate:');
    fprintf(' %.2f', convergence_rate);
    fprintf('\n');
end
