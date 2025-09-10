
function F = g(xt, p) % xt, each rows correspond each time, columns are different dofs

    % [N, M] = size(xt);
    % F.F = zeros(size(xt));
    % FF = zeros(size(xt))';
    % fc = p.fc;

    % Z1 = zeros(1, M)';
    % Z2 = zeros(2, 4);
    
    for j = 1:2
        for i = 1:256
            % x = xt(i, :);
            % Fi = gf(x, p);
            Fi = gf(xt(1, :), p);
            % F.F(i, :) = Fi.F';
            % F.F(i, :) = Z1;
            % FF(:, i) = Z1;
            % p.w = Fi.w; % update w
            % p.w = Z2;
            % w = Z2;
        end
    end
    % F.w = Fi.w; % update w to outside
    % w = Z2;
    F = 0;
end