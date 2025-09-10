clear
load("fc_column.mat");
% load("fc_row.mat");

N = 2^8;
Nx = 4;
% xtest = rand(N, 3 * Nx);
xtest = rand(N, 3 * Nx)';

tic;
for i = 1:100000
    % [F, w] = gf(xtest(1,:), fc);
    [F, w] = gf(xtest(:,1), fc);
end
toc;

