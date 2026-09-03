alpha_funs = {
    @(x) 0.9 - 0.1*sin(abs(pi*x))        
};


N_values = [16, 32, 64, 128, 256, 512];

[error_max, convergence] = error_analysis_MMS(N_values,alpha_fun{1});


