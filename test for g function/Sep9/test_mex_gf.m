clear
load("fc.mat");
N = 2^8;
Nx = 4;
xtest = rand(N, 3 * Nx);

% profile on
tic;
for i = 1:100000
    F = gf_mex(xtest(1,:), fc);
end
toc;
% profile off
% p = profile('info');
% profview(0,p);