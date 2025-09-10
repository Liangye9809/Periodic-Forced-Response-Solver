

function [F, w] = g(xt, fc) % xt, each rows correspond each time, columns are different dofs

    [N, ~] = size(xt);
    x = xt';
    F = zeros(size(x));
    for j = 1:2
        for i = 1:N
            [F(:, i), w] = gf(x(:, i), fc);
            fc.w = w;
        end
    end
    
end