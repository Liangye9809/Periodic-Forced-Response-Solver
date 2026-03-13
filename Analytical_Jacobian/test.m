skip = [1 2; 2 1; 3 1; 3 2];
for i = 1:3
    for j = 1:3
        if ismember([i j], skip, 'rows')
            continue;
        end
        disp([i,j]);
    end
end
%% plot 2D tangential relative displacement

T = 1;
N = 256;
phy = pi / 2;
dt = T / N;
t = 0:dt:10*pi;

figure;
plot(sin(t),cos(t), 'k-'), hold on;
grid on;
for Uy0 = 0.1:0.2:1
% Uy0 = 1;
    Ux0 = 2.5 * Uy0;
    
    Ux = Ux0 * sin(2*pi*t);
    Uy = Uy0 * sin(2*pi*t + phy);
    % figure;
    plot(Ux, Uy), hold on;

% drawnow;
% pause(0.05);
end
xlim([-3,3]);
ylim([-3,3]);
pbaspect([1 1 1]);

%% find transition time, from slip or separation to stick
a = [1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 -1 -1 0 1  1  1]'; % 0 is slip or separation. 1 is stick
% a(:, 1) = Mft_A(1,1,:);
% b(:, 1) = Mft_A(2,1,:);
% c(:, 1) = Mft_A(3,1,:);
% ll = 1:size(a,1);
N = length(a);
record = zeros(size(a));
p_plus = 0;
keep = 0;
for i = 1:N

    if i == 1 && a(1) ~=0 
        for j = 1:N
            k = N - j + 1;
            if a(k) == 0
                p_plus = mod(k, N) + 1;
                keep = 1;
                break;
            end
        end
    elseif a(i) ~= 0 && keep == 0
        p_plus = i;
        keep = 1;
    end

    if a(i) == 0 && keep == 1
        keep = 0;
        p_plus = 0;
    end

    record(i) = p_plus;


end
% c_xn = a .* b(mod(record-2, N) + 1);
% 
% conditi = [a, b, c_xn, record];
conditi = [a, record];

%% find transition time, from slip or separation to stick, in the middle
a = [0 1 1 0 0 0 0 0 1 1 1 1 1 1 0 -1 -1 0 1  0  0]'; % 0 is slip or separation. 1 is stick
% a(:, 1) = Mft_A(1,1,:);
% b(:, 1) = Mft_A(2,1,:);
% c(:, 1) = Mft_A(3,1,:);
ll = 1:size(a,1);
N = length(a);
record = zeros(size(a));
record_plus = zeros(size(a));
p_plus = 1; % set default value to 1 not 0, to avoid index of vector irreasonable in code
p_minus = 1;
p = 1;
keep = 0;
for i = 1:N

    if i == 1 && a(1) ~=0 % if begin with stick, means transition time located previously 
        for j = 1:N % backwards loop to find slip to stick transition
            k = N - j + 1;
            if a(k) == 0
                p_plus = mod(k, N + 1) + 1; % t+ transition instant, stick. add last point to create close the cycle [1, N+1]
                p_minus = mod(k - 1, N) + 1; % t- transition instant, slip
                p = 0.5 * (p_plus + p_minus); % middle transition instant.
                keep = 1;
                break;
            end
        end
    elseif a(i) ~= 0 && keep == 0
        p_plus = i; % t+ transition instant.
        p_minus = i - 1; % t- transition instant.
        p = 0.5 * (p_plus + p_minus); % middle transition instant.
        keep = 1;
    end

    if a(i) == 0 && keep == 1
        keep = 0;
        p_plus = 1;
        p_minus = 1;
        p = 1;
    end

    record_plus(i) = p_plus;
    record(i) = p_minus;


end
dxt = [1:21]';
dxt_p = dxt(record);
conditi = [a, record, record_plus, dxt, dxt_p];

%% test dXinFourier
N = 1024;
H = 200;
t = (0:(N-1)) * 2 * pi / N;
t = t';
% x = t .^ 3;
% dx = 3 * t .^ 2;
x = sin(sin(t));
dx = cos(sin(t)) .* cos(t);
[E EH] = fft_matrices(N, H);
X = EH * x;
dX = dXinFourier(X, H);
dxR = E * dX;

figure;
plot(t, dx, 'b-'), hold on;
plot(t, dxR, 'r--'), hold on;
grid on;

%%
dxt = [1:32]';
for i = 1:32
    disp(dxt(mod(i - 2, 32) + 1));
end

%%
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
H = 5;
Nx = 16;
n = 3 * Nx * (2 * H + 1); % size of matrix
[E, EH] = fft_matrices(N, H);
JNLt_rand = rand(N, 2 * H + 1);
Loop = 100;

n_range = 4; % stick-slip-stick-slip
% n_W = 3; % number of W, dt/dxt, dt/dxn, dfn/dxn
tic;
for i = 1:Loop * (32 * H) * n_range
    a = cos(rand(1));
end
T_W_cos = toc

n_d = 3; % number of dF, dt/dxt, dt/dxn, dfn/dxn
tic;
for i = 1:Loop * n_d
    % a = EH * JNLt_rand;
    a = EH * rand(N, 2 * H + 1);
    b = E .* rand(N, 1);
    % for j = 1:(2 * H * N)
    %     b = cos(rand(1));
    % end
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

%% analytcial and numerical Fourier comparison
t1 = 1;
t2 = 2;
H = 5;
N = 128;
epsW = [];
epsw = [];
for iN = 1:10
    N = 2^(iN + 7);
    for iH = 1:20
        H = iH;
        [E, EH] = fft_matrices_t1_t2(N, H, t1, t2);
        [W_A, w_a] = fW(t1, t2, H);
        W_N = EH * E;
        w_n = EH * ones(N, 1);
        epsW(iN, iH) = norm(W_A(2:end, :) - W_N(2:end, :)) / norm(W_A(2:end, :));
        epsw(iN, iH) = norm(w_a - w_n) / norm(w_a);

    end
end

%%
a = zeros(32, 2);
a([10, 20], 1) = 1;
a([10, 20], 2) = 3.14;
b = find(a(:,1) > 0)
c = a(b, 2)