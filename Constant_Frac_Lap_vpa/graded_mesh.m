function x = graded_mesh(a, b, N, alp)

    a   = vpa(a);
    b   = vpa(b);
    alp = vpa(alp);

    two = vpa(2);
    m   = vpa(N)/two;

    r = two*(two-alp)/alp;

    j_double = 0:N;
    j = vpa(j_double);

    x = sym(zeros(1, N+1));

    left_indices = (j_double <= N/2);

    x(left_indices) = a + (b-a)/two .* ...
        (j(left_indices)/m).^r;

    x(~left_indices) = b - (b-a)/two .* ...
        (two-j(~left_indices)/m).^r;
end
