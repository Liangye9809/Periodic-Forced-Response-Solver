function [W, w] = fW1_sym(t1, t2, H)

    if t1 > t2
        t2 = t2 + 2*pi;
    end

    W = zeros(2 * H + 1, 2 * H + 1);
    w = zeros(2 * H + 1, 1);

    W(1, 1) = t2 - t1;

    for i = 1:H
        % First column / row (unchanged)
        W(2 * i, 1) = (1 / i) * (sin(i * t2) - sin(i * t1));
        W(2 * i + 1, 1) = - (1 / i) * (cos(i * t2) - cos(i * t1));

        W(1, 2 * i) = 2 * W(2 * i, 1);
        W(1, 2 * i + 1) = 2 * W(2 * i + 1, 1);

        % Diagonal block (i == j)
        a = (1 / (2*i)) * (cos((2*i) * t2) - cos((2*i) * t1));
        c = (1 / (2*i)) * (sin((2*i) * t2) - sin((2*i) * t1));

        W(2 * i, 2 * i) = c + (t2 - t1);
        W(2 * i + 1, 2 * i) = - a;
        W(2 * i, 2 * i + 1) = - a;
        W(2 * i + 1, 2 * i + 1) = - c + (t2 - t1);

        % Off-diagonal (j < i only)
        for j = 1:i-1
            a = (1 / (i + j)) * (cos((i + j) * t2) - cos((i + j) * t1));
            c = (1 / (i + j)) * (sin((i + j) * t2) - sin((i + j) * t1));

            b = (1 / (i - j)) * (cos((i - j) * t2) - cos((i - j) * t1));
            d = (1 / (i - j)) * (sin((i - j) * t2) - sin((i - j) * t1));

            % Compute once
            W(2 * i, 2 * j) = c + d;
            W(2 * i + 1, 2 * j) = - a - b;
            W(2 * i, 2 * j + 1) = - a + b;
            W(2 * i + 1, 2 * j + 1) = - c + d;

            % Symmetry copy (same as fW2 idea)
            W(2 * j, 2 * i) = W(2 * i, 2 * j);
            W(2 * j, 2 * i + 1) = W(2 * i + 1, 2 * j);
            W(2 * j + 1, 2 * i) = W(2 * i, 2 * j + 1);
            W(2 * j + 1, 2 * i + 1) = W(2 * i + 1, 2 * j + 1);
        end
    end

    W = (1 / (2 * pi)) .* W;
    w = 2 .* W(:, 1);

end