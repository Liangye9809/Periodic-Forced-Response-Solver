function [J, Fi, xi] = finite_diff_jac_F(f, x)
    
    [f_i, x_i] = f(x);
    
    n = length(x);   
    m = length(f_i);   
    J = zeros(m, n);
    h = 1e-8;
    
    Fi = f_i;
    xi = x_i;
    for j = 1:n
        x_ip = x;
        x_ip(j) = x_ip(j) + h;
        [f_ip, x_t] = f(x_ip);
        Fi = [Fi, f_ip];
        xi = [xi, x_t];
        J(:, j) = (f_ip - f_i) / h;
    end
end
