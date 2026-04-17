function c_row = c_vector(tau, H)
    c_row = zeros(1, 2 * H + 1);
    c_row(1) = 0.5;
    % c_row(1) = 1;
    for i = 1:H
        c_row(2 * i) = cos(i * tau);
        c_row(2 * i + 1) = sin(i * tau);
    end

    
end

