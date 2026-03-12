function [E, EH] = fft_matrices_t1_t2(N, H, t1, t2)
            E = zeros(N, 2 * H + 1);
            E(:, 1) = 1 / 2;

            if t1 > t2
                t2 = t2 + 2 * pi; % make sure t2 > t1
            end

            t = (0:(N - 1))';
            t = t1 + ((t2 - t1) / (N - 1)) * t;
            for k = 1:H
                E(:, 2 * k) = cos(t * k);
                E(:, 2 * k + 1) = sin(t * k);
            end
            
            EH = E';
            EH(1,:) = 1;
            
            EH = EH * (2 / N) * (t2 - t1) / (2 * pi);
end
