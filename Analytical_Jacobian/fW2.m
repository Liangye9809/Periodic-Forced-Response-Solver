function [W, w] = fW2(t1, t2, H)
    if t1 > t2
        t2 = t2 + 2*pi;
    end
    W = zeros(2 * H + 1, 2 * H + 1);
    w = zeros(2 * H + 1, 1);

    

    cost = zeros(2 * H, 1);
    sint = zeros(2 * H, 1);

    for i = 1:2*H
        cost(i) = (cos(i * t2) - cos(i * t1)) / i;
        sint(i) = (sin(i * t2) - sin(i * t1)) / i;
    end

    % i = [1:2*H]';
    % cost = (cos(i * t2) - cos(i * t1)) ./ i;
    % sint = (sin(i * t2) - sin(i * t1)) ./ i;

    for i = 1:H
        W(2 * i, 1) = sint(i);
        W(2 * i + 1, 1) = - cost(i);
        W(1, 2 * i) = 2 * W(2 * i, 1);
        W(1, 2 * i + 1) = 2 * W(2 * i + 1, 1);

        W(2 * i, 2 * i) = sint(2 * i) + (t2 - t1);
        W(2 * i + 1, 2 * i) = - cost(2 * i);
        W(2 * i, 2 * i + 1) = - cost(2 * i);
        W(2 * i + 1, 2 * i + 1) = - sint(2 * i) + (t2 - t1);
        for j = 1:i-1
            W(2 * i, 2 * j) = sint(i + j) + sint(i - j);
            W(2 * i + 1, 2 * j) = - cost(i + j) - cost(i - j);
            W(2 * i, 2 * j + 1) = - cost(i + j) + cost(i - j);
            W(2 * i + 1, 2 * j + 1) = - sint(i + j) + sint(i - j);

            W(2 * j, 2 * i) = W(2 * i, 2 * j);
            W(2 * j, 2 * i + 1) = W(2 * i + 1, 2 * j);
            W(2 * j + 1, 2 * i) = W(2 * i, 2 * j + 1);
            W(2 * j + 1, 2 * i + 1) = W(2 * i + 1, 2 * j + 1);
        end
    end
    
    W(1, 1) = t2 - t1;
    
    % for i = 1:H
    %     W(2 * i, 1) = sint(i);
    %     W(2 * i + 1, 1) = - cost(i);
    %     W(1, 2 * i) = 2 * W(2 * i, 1);
    %     W(1, 2 * i + 1) = 2 * W(2 * i + 1, 1);
    % 
    %     W(2 * i, 2 * i) = sint(2 * i) + (t2 - t1);
    %     W(2 * i + 1, 2 * i) = - cost(2 * i);
    %     W(2 * i, 2 * i + 1) = - cost(2 * i);
    %     W(2 * i + 1, 2 * i + 1) = - sint(2 * i) + (t2 - t1);
    % end


    % for i = 1:H
    %     for j = 1:i
    %         if i == j
    %             W(2 * i, 2 * i) = sint(2 * i) + (t2 - t1);
    %             W(2 * i + 1, 2 * i) = - cost(2 * i);
    %             W(2 * i, 2 * i + 1) = - cost(2 * i);
    %             W(2 * i + 1, 2 * i + 1) = - sint(2 * i) + (t2 - t1);
    %         else
    %             W(2 * i, 2 * j) = sint(i + j) + sint(i - j);
    %             W(2 * i + 1, 2 * j) = - cost(i + j) - cost(i - j);
    %             W(2 * i, 2 * j + 1) = - cost(i + j) + cost(i - j);
    %             W(2 * i + 1, 2 * j + 1) = - sint(i + j) + sint(i - j);
    %         end
    %     end
    % end

    W = (1 / (2 * pi)) .* W;

    w = 2 .* W(:, 1);
    

end
