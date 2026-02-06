function [J, w, Mft] = finite_diff_jac(f, x, h, order)
    
    % f_i = f(x);
    [f_i, w, Mft] = f(x);
    
    n = length(x);   
    m = length(f_i);   
    J = zeros(m, n);
    % h = 1e-8;
    % h = 1e-6;
    if order == 1
        for j = 1:n
            x_ip = x;
            x_ip(j) = x_ip(j) + h;
            % f_ip = f(x_ip);
            [f_ip, w, Mft] = f(x_ip);
            J(:, j) = (f_ip - f_i) / h;
        end
    elseif order == 2
        for j = 1:n
            x_ip = x;
            x_ip(j) = x_ip(j) + h;
            % f_ip = f(x_ip);
            [f_ip, w, Mft] = f(x_ip);

            x_im = x;
            x_im(j) = x_im(j) - h;
            % f_im = f(x_im);
            [f_im, w, Mft] = f(x_im);
            J(:, j) = (f_ip - f_im) / (2 * h);
        end
    end
end
