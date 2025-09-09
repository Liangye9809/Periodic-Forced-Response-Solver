% x(t) are in time domain in conlumn vector, 
% number of columns is number of DOF, N*Nc matrix
% only in contact part

% function F = g(xt, p) 
%     beta = p.Nondimention.beta;
%     alpha = p.Nondimention.alpha;
%     omega02 = p.Nondimention.omega02;
%     F.F = (beta^4 / (alpha^2*omega02)) * 5e16 * xt.^3; % 2.79e10 * xt^3
%     % F.F = (beta^4 / (alpha^2*omega02)) * xt.^3;
% 
%     % F.F = 5e16 .* xt.^3;
% 
% 
%     % F.F = 0*xt;
%     % F.F = 1 .* xt.^3;
%     F.w = p.fc.w;
% end

%% for mex function
% function F = g(xt, p) % xt, each rows correspond each time, columns are different dofs
% 
%     [N, ~] = size(xt);
%     F.F = zeros(size(xt));
%     % convert fc to struct type
%     fc.kn = p.kn;
%     fc.xn0 = p.xn0;
%     fc.mu = p.mu;
%     fc.kt = p.kt;
%     fc.w = p.w;
% 
%     for j = 1:2
%         for i = 1:N
%             x = xt(i, :);
%             [FF, fc.w] = gf_mex(x, fc);
%             F.F(i, :) = FF;
%         end
%     end
%     F.w = fc.w;
% 
% end

%% matlab only
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