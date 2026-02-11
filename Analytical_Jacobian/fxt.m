function x = fxt(t)
    
    % x(:, 1) = exp(sin(t));
    x(:, 1) = exp(cos(t));
    x(:, 2) = 1 ./ (1 + (cos(t)).^2);
    x(:, 3) = sin(sin(t));

end