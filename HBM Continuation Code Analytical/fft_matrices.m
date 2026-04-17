function [E, EH] = fft_matrices(N, H)
    E = zeros(N, 2 * H + 1);
    E(:, 1) = 1 / 2;

    t = (0:(N - 1))';
    
    for h = 1:H
        E(:, 2 * h) = cos(2*pi/N * t * h);
        E(:, 2 * h + 1) = sin(2*pi/N * t * h);
    end
    
    EH = E';
    EH(1,:) = 1;
    
    EH = EH * (2/N);
end