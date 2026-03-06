function [W, w] = fW(t1, t2, H)
    W = zeors(2 * H + 1, 2 * H + 1);
    w = zeros(2 * H + 1, 1);
    W(1, 1) = t2 - t1;
    for i = 1:H
        W(1, 2 * i) = 1 / i * (sin(i * t2) - sin(i * t1));
        W(1, 2 * i + 1) = - 1 / i *(cos(i * t2) - cos(i * t1));
        W(2 * i, 1) = 2 * W(1, 2 * i);
        W(2 * i + 1, 1) = 2 * W(1, 2 * i + 1);

        for j = 1:H
            a = 1 / (i + j) * (cos((i + j) * t2) - cos((i + j) * t1));
            c = 1 / (i + j) * (sin((i + j) * t2) - sin((i + j) * t1));
            if i == j
                W(2 * i, 2 * i) = c + (t2 - t1);
                W(2 * i + 1, 2 * i) = - a;
                W(2 * i, 2 * i + 1) = - a;
                W(2 * i + 1, 2 * i + 1) = - c + (t2 - t1);
            else
                b = 1 / (i - j) * (cos((i - j) * t2) - cos((i - j) * t1));
                d = 1 / (i - j) * (sin((i - j) * t2) - sin((i - j) * t1)); 
                W(2 * i, 2 * i) = c + d;
                W(2 * i + 1, 2 * i) = - a - b;
                W(2 * i, 2 * i + 1) = - a + b;
                W(2 * i + 1, 2 * i + 1) = - c + d;
            end
        end
    end
    W = (1 / (2 * pi)) .* W;

    w = W(:, 1);

end

