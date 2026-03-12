clear

% a = [0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0]';
% b = [-1,-1,-1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,-1,-1]';
% c = [1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1]';

a = 0 .* [1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0]';
b = 0 .* [0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0]';
c = 0 .* [1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0]';

Mf = [a, b, c];
N = size(a, 1);
[C_ss, C_gs, C_sp, C_sm, C_c] = get_integral_time_position(a, b, c);
n_m = mod(ceil(C_gs(1, 1)) - 1, N) + 1;
n_p = mod(ceil(C_gs(1, 1) + 1) - 1, N) + 1;
%%
i = 1;
while i < 10
    disp(i)
    i = i + 1;
end

c_row = c_vector(0.5 * pi, 5);

%%
clear
clc
A = rand(5, 6)
rankA = rank(A)
[U, S, V] = svds([A; zeros(1, 6)], 1, 'smallest');
A_times_V = A * V

% U_times_A = U' * A

nullA = null(A)

A_times_nullA = A * nullA

%%
clear
clc
N = 256;
H = 10;
Nx = 16;

h = 10^(-4);
order = 1;
dt = 2 * pi / N;
t = (0:(N-1)) * 2 * pi / N;
t = t';
% xn = ones(N, 1);
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xn = 2 * exp(cos(t)) - 0.75; % separation to slip
xt = 2 * sin(sin(t)); % slip to stick
% xt = sin(sin(t)); % pure stick
x = [xt, xn];



kt = 1;
kn = 2;
mu = 0.5;
w =  -1;
xn0 = 0; % normal pre-displacement

nloop = 2;

[E, EH] = fft_matrices(N, H);
X = EH * x;
dX = dXinFourier(X, H);
dx = E * dX;
X = [EH * x(:, 1), X];
X = X(:);
[Ft, wt, Mft] = gf_2dofs(x, kn, xn0, mu, kt, w, nloop);

Loop = 10000;

tic;
for i = 1:Loop
    [JNL_A, ~] = HBMJACOB_analytical_gf_2dofs_2(dx, kn, mu, kt, H, N, Mft);
end
T_A = toc

tic;
for i = Loop
    JNL_W = HBMJACOB_analytical_W_gf_2dofs_2(dx, kn, mu, kt, H, N, Mft);
end
T_W = toc

T_A / T_W
function dX = dXinFourier(X, H)
    dX = zeros(size(X));
    for i = 1:H
        dX(2 * i, :) =  i .* X(2 * i + 1, :);
        dX(2 * i + 1, :) =  -i .* X(2 * i, :);
    end

end

%%
clear
% clc
N = 256;
H = 10;
Nx = 16;
n = Nx * (2 * H + 1); % size of matrix
[E, EH] = fft_matrices(N, H);
JNLt_rand = rand(N, 2 * H + 1);
Loop = 100;

n_range = 4; % stick-slip-stick-slip
% n_W = 3; % number of W, dt/dxt, dt/dxn, dfn/dxn
tic;
for i = 1:Loop * (8 * H * H) * n_range
    a = cos(rand(1));
end
T_W_cos = toc

n_d = 3; % number of dF, dt/dxt, dt/dxn, dfn/dxn
tic;
for i = 1:Loop * n_d
    % a = EH * JNLt_rand;
    a = EH * rand(N, 2 * H + 1);
    for j = 1:(2 * H * N)
        b = cos(rand(1));
    end
end
T_A_EH = toc
%% null
clear
clc
N = 256;
H = 10;
Nx = 16;
n = Nx * (2 * H + 1); % size of matrix
A = rand(n, n + 1);
Fx = A(:, 1:n);
Fl = A(:, end);
Loop = 1000;


tic;
for i = 1:Loop
    a = null(A);
end
T_null = toc


tic;
for i = 1:Loop
    a = svds([A; zeros(1, n + 1)], 1, 'smallest');
end
T_svds = toc

tic;
for i = 1:Loop
    b = - (Fx\Fl);
    if norm(b) > 1e10
        a = null(A);
    else
        a = [b; 1];
        a = a / norm(a);
    end
end
T_direc = toc

