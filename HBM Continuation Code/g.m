% x(t) are in time domain in conlumn vector, 
% number of columns is number of DOF, N*Nc matrix
% only in contact part

% function F = g(xt, p) 
%     F = 5e16 .* xt.^3;
%     % F = sin(x);
%     % F = 0*xt;
%     % F = xt.^3;
% end

function F = g(xt, p) % xt, each rows correspond each time, columns are different dofs

    [N, ~] = size(xt);
    F = zeros(size(xt));
    % pc = p.fc;

    for j = 1:2
        for i = 1:N
            x = xt(i, :);
            Fi = gf(x, p);
            F(i, :) = Fi';
        end
    end


end