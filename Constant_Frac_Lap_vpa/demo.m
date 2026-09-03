clear;
clc;


digits(80);

N_values = [32,64,128,256,512];

alp_values = sym('0.3');

[error_max, convergence] = ...
    error_analysis(N_values, alp_values);
