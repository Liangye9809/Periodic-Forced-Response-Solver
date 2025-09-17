function J = finite_diff_jac(f, x)
    
    f_i = f(x);
    
    n = length(x);   
    m = length(f_i);   
    J = zeros(m, n);
    h = 1e-8;
    
    for j = 1:n
        x_ip = x;
        x_ip(j) = x_ip(j) + h;
        f_ip = f(x_ip);
        J(:, j) = (f_ip - f_i) / h;
    end
end
