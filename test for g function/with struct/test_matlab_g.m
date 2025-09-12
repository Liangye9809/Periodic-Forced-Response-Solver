clear
load("fc_column.mat");
N = 2^8;
Nx = 4;
xtest = rand(N, 3 * Nx);

profile on
tic;
for i = 1:1000
    [F, w] = g(xtest, fc);
end
toc;
profile off
p = profile('info');
profview(0,p);
