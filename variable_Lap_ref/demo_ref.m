

alpha_funs = {
    @(x) 0.9 - 0.3*abs(x)        
};


N_values = [16, 32, 64, 128, 256, 512];


[error_max, convergence] = error_analysis_ref(N_values, alpha_funs{1});