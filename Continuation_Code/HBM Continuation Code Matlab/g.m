
function [F, w] = g(xt, kn, xn0, mu, kt, w_in) % xt, each rows correspond each time, columns are different dofs

    [N, ~] = size(xt);
    x = xt';
    F = zeros(size(x));

    for j = 1:2
        for i = 1:N
            [F(:, i), w] = gf(x(:, i), kn, xn0, mu, kt, w_in);
            w_in = w; 
        end
    end
    
    F = F';
end