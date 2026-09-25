%% Coulomb friction of dummy fucntion 2 dofs
clear
clc
% close all
eps = [];
h = 10^(-7);
order = 1;
h_con = [];
N = 16;
H = 3;
dt = 2 * pi / N;
t = (0:(N-1)) * 2 * pi / N;
t = t';
% xn = ones(N, 1);
% xn = - 4 * sin(sin(t)) + 1; % separation to stick
% xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
% xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xn = 2 * exp(cos(t)) - 0.75; % separation to slip
% xt = 2 * sin(sin(t)); % slip to stick

% xt = 1.05  * sin(2 .* exp(cos(t))); % tangent case1
% xt = 1.00 * sin(sin(t)) ./ sin(1); % tangent case2
% xt = 0.5 * sin(sin(t)) ./ sin(1) + 0.5; % tangent case3 only one side

% plot w for the paper (kt = 1, kn = 2, mu = 0.5)
% pure stick
% xn = 2*ones(N, 1); % pure stick
% xt = sin(sin(t)); % pure stick
% simple x
% xn = 2*ones(N, 1); % pure stick
% xt = sin(t); % pure stick

% slip to stick
% xn = 2 * exp(cos(t)) - 0.5; % slip to stick
% xt = 2 * sin(sin(t)); % slip to stick
% simple x
% xn = 2.5 * cos(t) + 3; % slip to stick
% xt = 3 * sin(t); % slip to stick

% gap to stick
xn = - 4 * sin(sin(t)) + 1; % separation to stick
xt = 2 * exp(cos(t + 1)) - 3; % separation to stick
% simple x
% xn = - 10 * cos(t) + 3; % separation to stick
% xt = - sin(t); % separation to stick
% xn = - 1 * cos(t) + 0.5; % separation to stick % (kt = 1, kn = 500, mu = 0.8)
% xt = 2 * cos(t); % separation to stick
x = [xt, xt, xn];

[E, EH] = fft_matrices(N, H);
X = EH * x;
xpr = E * X;
dX = dXinFourier(X, H);
dx = E * dX;
dxt1 = dx(:,1);
dxn = dx(:,3);

kt = [1;1];
kn = 2;
mu = [0.5;0.5];
w =  [0;0];
xn0 = 0; % normal pre-displacement

nloop = 2;
[Fti, wi, flag] = g(x, kn, xn0, mu, kt, w, nloop);

segments = get_integral_time_position(flag(1,1,end - N + 1:end), xt, xn, dxt1, dxn, kt(1), kn, mu(1), Fti(end - N + 1:end, 1), H)
S{1} = segments;
S{2} = segments;
JNL = JNL_Analytical(S, H, kt, kn, mu);

function dX = dXinFourier(X, H)
    dX = zeros(size(X));
    for i = 1:H
        dX(2 * i, :) =  i .* X(2 * i + 1, :);
        dX(2 * i + 1, :) =  -i .* X(2 * i, :);
    end

end

%% compare the F and Jacobian
clear
clc
close all
GetDataTest
dxct = get_dxt(Xc, E, Nx);

[Fti, wi, flagi] = g(xct, kn, xn0, mu, kt, w, nloop, dxct); 
Ft = Fti(end - N + 1:end, :); % the last period
flag = flagi(:, :, end - N + 1:end); % the last period
S = get_all_segments(flag, xct, dxct, kt, kn, mu, Ft, H);

J_a = JNL_Analytical(S, H, kt, kn, mu);
F_a = get_Analytical_F_Fourier(S, gxp, kt, kn, mu, H, Xc);
J_a_pre = JNL_Analytical_pre(xct, flag, H, N, kt, kn, mu);

dft = FFT_improve(Ft, xct, flag, kt, kn, mu);
hndn_pre = EH * (Ft - gxp');
F_N_pre = hndn_pre(:);
hndn = EH * (Ft - gxp' + dft);
F_N = hndn(:);

J_N = finite_diff_jac(@(Xc) fftF(Xc, H, Nx, kn, xn0, mu, kt, w, nloop, gxp, E, EH, N, dxct), Xc);



function F_Fourier = fftF(Xc, H, Nx, kn, xn0, mu, kt, w, nloop, gxp, E, EH, N, dxct)
    XcM = zeros(2 * H + 1, 3 * Nx);
    for i = 1:3*Nx
        XcM(:, i) = Xc((2 * H + 1) * (i - 1) + 1:(2 * H + 1) * i);
    end
    xct = E * XcM;
    [Fti, wi, flag] = g(xct, kn, xn0, mu, kt, w, nloop, dxct); 
    flag_ = flag(:, :, end - N + 1:end);
    if sum(ismember([-1, 1, 0], flag_)) > 0 % slip gap appear
        Ft_in = Fti(end - N + 1:end, :);

        dft = FFT_improve(Ft_in, xct, flag_, kt, kn, mu);

        Fti(end - N + 1:end, :) = Fti(end - N + 1:end, :) + dft;
        
    end
    Ft = Fti(end - N + 1:end, :) - gxp';
    F_Fourier = EH * Ft;
    F_Fourier = F_Fourier(:);
end

epsF_pre = norm(F_N_pre - F_a) / norm(F_a)
epsF = norm(F_N - F_a) / norm(F_a)

epsF_pre1 = norm(F_N_pre(1:(2 * H + 1)) - F_a(1:(2 * H + 1))) / norm(F_a(1:(2 * H + 1)))
epsF1 = norm(F_N(1:(2 * H + 1)) - F_a(1:(2 * H + 1))) / norm(F_a(1:(2 * H + 1)))

epsF_pre2 = norm(F_N_pre((2 * H + 2):(4 * H + 2)) - F_a((2 * H + 2):(4 * H + 2))) / norm(F_a((2 * H + 2):(4 * H + 2)))
epsF2 = norm(F_N((2 * H + 2):(4 * H + 2)) - F_a((2 * H + 2):(4 * H + 2))) / norm(F_a((2 * H + 2):(4 * H + 2)))

epsF_pre3 = norm(F_N_pre((4 * H + 3):end) - F_a((4 * H + 3):end)) / norm(F_a((4 * H + 3):end))
epsF3 = norm(F_N((4 * H + 3):end) - F_a((4 * H + 3):end)) / norm(F_a((4 * H + 3):end))

epsJ_analytical_numerical = norm(J_N - J_a) / norm(J_N)
epsJ_analytical_pre_new = norm(J_a - J_a_pre) / norm(J_a)
epsJ_analytical_pre_numerical = norm(J_N - J_a_pre) / norm(J_N)

flagT1(:, 1) = flag(1, 1, :);
S1_pre = get_integral_time_position_pre(flagT1);
flagT2(:, 1) = flag(2, 1, :);
S2_pre = get_integral_time_position_pre(flagT2);

figure % forces
plot(mu(1) * (Ft(:, 3)), 'LineWidth', 2, 'LineStyle','-', 'DisplayName', '$\mu Fn$'), hold on;
plot(-mu(1) * (Ft(:, 3)), 'LineWidth', 2, 'LineStyle','-', 'DisplayName', '$-\mu Fn$'), grid on;
plot(Ft(:, 1), 'LineWidth', 2, 'LineStyle','-', 'DisplayName', '$Ft$'), hold on;
legend('show')

figure
plot((xct(:, 1)) - mu(1) * (Ft(:, 3)) / kt(1), 'LineWidth', 2, 'LineStyle','-', 'DisplayName', '$w^-$'), hold on;
plot((xct(:, 1)) + mu(1) * (Ft(:, 3)) / kt(1), 'LineWidth', 2, 'LineStyle','-', 'DisplayName', '$w^+$'), grid on;
wplot(:, 1) = wi(1,1,end - N +1:end);
plot(wplot, 'LineWidth', 2, 'LineStyle','-', 'DisplayName', '$w$'), hold on;
legend('show')

figure
plot(xct(:, 1), 'LineWidth', 2, 'LineStyle','-', 'DisplayName', '$xt1$'), hold on;
plot(xct(:, 3), 'LineWidth', 2, 'LineStyle','-', 'DisplayName', '$xn$'), hold on;
legend('show');
function JNL = JNL_Analytical_pre(x, flag, H, N, kt, kn, mu) % x is the size of N*3Nx

    Nx = size(flag, 2);
    JNL = zeros(3 * Nx * (2 * H + 1), 3 * Nx * (2 * H + 1));

    for i = 1:Nx
        flagT1(:, 1) = flag(1, i, :);
        flagT2(:, 1) = flag(2, i, :);
        JNLi = JNL_one_Nx(x(:, 3 * i - 2:3 * i), flagT1, flagT2, H, N, kt(:, i), kn(i), mu(:, i));
        indx1 = 3 * (i - 1) * (2 * H + 1) + 1;
        indx2 = 3 * i * (2 * H + 1);
        JNL(indx1:indx2, indx1:indx2) = JNLi;
    end
            

end


function JNLi = JNL_one_Nx(x, flagT1, flagT2, H, N, kt, kn, mu)

    JNLi = zeros(3 * (2 * H + 1), 3 * (2 * H + 1));

    segmentsT1 = get_integral_time_position_pre(flagT1);
    [dF1dX1, dF1dXn, dFndXn] = get_dFdX(segmentsT1, H, N, kt(1), kn, mu(1), x(:, 1), x(:, 3));
    
    segmentsT2 = get_integral_time_position_pre(flagT2);
    [dF2dX2, dF2dXn, ~] = get_dFdX(segmentsT2, H, N, kt(2), kn, mu(2), x(:, 2), x(:, 3));

    JNLi(1:2 * H + 1, 1:2 * H + 1) = dF1dX1;
    JNLi(1:2 * H + 1, 2 * (2 * H + 1) + 1:end) = dF1dXn;

    JNLi((2 * H + 1) + 1:2 * (2 * H + 1), (2 * H + 1) + 1:2 * (2 * H + 1)) = dF2dX2;
    JNLi((2 * H + 1) + 1:2 * (2 * H + 1), 2 * (2 * H + 1) + 1:end) = dF2dXn;

    JNLi(2 * (2 * H + 1) + 1:end, 2 * (2 * H + 1) + 1:end) = dFndXn;
end


function segments = get_integral_time_position_pre(flag)
    % FIND_CYCLIC_SEGMENTS Find cyclic time ranges for each segment in flag vector
    % Time axis: t(i) = (i-1)/N * 2pi, i = 1..N
    % Boundaries are at midpoints between transitions
    
    flag = flag(:);
    N = length(flag);
    
    % --- Find transition points (cyclic) ---
    % diff detects where value changes; also check wrap-around
    diffs = [diff(flag); flag(1) - flag(end)];  % length N, last element is wrap
    trans_idx = find(diffs ~= 0);               % indices WHERE change happens
    % Boundary midpoint index (fractional): between trans_idx and trans_idx+1 (cyclic)
    n_trans = length(trans_idx);
    
    if n_trans == 0
        % fprintf('Signal is constant: value = %d over [0, 2pi]\n', flag(1));
        segments = struct('value', flag(1), 't_start', 0, 't_end', 2*pi, 'index_start', 1);
        return;
    end
    
    % Boundary midpoints (fractional 1-based index, cyclic)
    boundary_idx = trans_idx + 0.5;  % midpoint between last-of-old and first-of-new
    
    % --- Build segments ---
    % Each segment starts at one boundary and ends at the next
    % Value of segment = flag at the first index AFTER the boundary
    
    segments = struct('value', {}, 't_start', {}, 't_end', {}, 'index_start', {});
    
    for k = 1:n_trans
        b_start = boundary_idx(k);                          % start boundary (index)
        b_end   = boundary_idx(mod(k, n_trans) + 1);        % next boundary (index, cyclic)
        
        % First sample index after b_start
        first_idx = mod(trans_idx(k), N) + 1;               % 1-based, cyclic
        val = flag(first_idx);
        
        % Convert boundary indices to time: t = (idx - 1) / N * 2pi
        t_start = (b_start - 1) / N * 2 * pi;
        t_end   = (b_end   - 1) / N * 2 * pi;
        
        % Handle cyclic wrap: if t_end <= t_start, it wraps around
        if t_end <= t_start
            t_end = t_end + 2 * pi;
        end
        
        segments(k).value       = val;
        segments(k).t_start     = t_start;
        segments(k).t_end       = t_end;
        segments(k).index_start = mod(b_start - 0.5, N) + 1;
    end
end

function [dFdX, dFdXn, dFndXn] = get_dFdX(segmentsT, H, N, kt, kn, mu, xt, xn)
    dFdX   = zeros(2 * H + 1, 2 * H + 1);
    dFdXn  = zeros(2 * H + 1, 2 * H + 1);
    dFndXn = zeros(2 * H + 1, 2 * H + 1);

    k = length(segmentsT);

    if k == 1 % pure stick or whole gap
        if segmentsT(1).value == 2 % pure stick
            dFdX   = kt .* eye(2 * H + 1);
            dFndXn = kn .* eye(2 * H + 1);
            return;
        elseif segmentsT(1).value == 0 % whole gap
            return;
        end
    end

    for i = 1:k
        t1 = segmentsT(i).t_start;
        t2 = segmentsT(i).t_end;
        switch segmentsT(i).value
            case 2 % stick after slip or gap
                [MW, Mw] = fW(t1, t2, H);
                c_vec = c_vector(t1, H);
                dFdX  = dFdX + kt .* MW - Mw * (kt .* c_vec);
                dFndXn = dFndXn + kn .* MW; % contact
                switch segmentsT(mod(i - 2, k) + 1).value

                    case 0 % gap to stick
                        i_p   = segmentsT(i).index_start; % stick start index
                        i_m   = mod(i_p - 2, N) + 1;  % previous index
                        dxdn  = (xt(i_p) - xt(i_m)) / (xn(i_p) - xn(i_m));
                        dFdXn = dFdXn + Mw * (kt .* dxdn .* c_vec);

                    case {-1, 1} % slip to stick
                        dFdXn = dFdXn + Mw * (segmentsT(mod(i - 2, k) + 1).value .* mu .* kn .* c_vec);
                end 

            case {-1, 1} % slip
                [MW, ~] = fW(t1, t2, H);
                dFdXn   = dFdXn + segmentsT(i).value .* mu .* kn .* MW;
                dFndXn  = dFndXn + kn .* MW; % contact
        end

    end

end


%%
for i = 1:29
    if ismember('data', h(i).DisplayName) 
        h(i).HandleVisibility = 'off';
    end
end

%%
clear
N = 32;
dt = 2 * pi / N;
t = (0:(N-1)) * 2 * pi / N;
t = t';
kn = 1;
kt = 0.5;
mu = 0.5;
w0 = 0;
xn = cos(t) + 1.5; 
xt = 2 * sin(t) + 0.2; 
[fn, ~, ft, ~, w] = TangentialForces_withPre(xn, mu, xt, kt, N, w0, kn);

N_ = 2^12;
dt_ = 2 * pi / N_;
t_ = (0:(N_-1)) * 2 * pi / N_;
t_ = t_';
xn_ = cos(t_) + 1.5; % slip to stick
xt_ = 2 * sin(t_) + 0.2; % slip to stick
[fn_, fnp, ft_, ftp, w_] = TangentialForces_withPre(xn_, mu, xt_, kt, N_, w0, kn);

% figure(100)
% plot(t, xn, 'bo-', LineWidth=2, DisplayName='$x_n$'), hold on, grid on;
% plot(t, fn, 'r*-', LineWidth=2, DisplayName='$f_n$'), hold on, grid on;
% legend('show')
% xticks([0:pi/2:2*pi]);
% xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
% set(gca,'TickLabelInterpreter','latex');
% xlim([0, 2*pi]);

figure(200)
plot(t_, mu * fnp, 'g-', LineWidth=2, DisplayName='$\mu f_n^{pre}$'), hold on, grid on;
plot(t_, -mu * fnp, 'g-', LineWidth=2, DisplayName='$\mu f_n^{pre}$', HandleVisibility='off'), hold on, grid on;
plot(t, mu * fn, 'bo-', LineWidth=2, DisplayName='$\mu f_n$'), hold on, grid on;
plot(t, -mu * fn, 'bo-', LineWidth=2, DisplayName='$\mu f_n$', HandleVisibility='off'), hold on, grid on;
plot(t_, ftp, 'r-', LineWidth=2, DisplayName='$f_t^{pre}$'), hold on, grid on;
plot(t, ft, 'k+-', LineWidth=2, DisplayName='$f_t$'), hold on, grid on;

legend('show')
xticks([0:pi/2:2*pi]);
xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
set(gca,'TickLabelInterpreter','latex');
xlim([0, 2*pi]);

w_m = xt_ - mu * fn_ / kt;
w_p = xt_ + mu * fn_ / kt;
figure(300)
plot(t_, w_m, 'r--', LineWidth=2, DisplayName='$w^{-}$'), hold on, grid on;
plot(t_, w_p, 'k--', LineWidth=2, DisplayName='$w^{+}$'), hold on, grid on;
plot(t, w, 'bo-', LineWidth=2, DisplayName='$w$'), hold on, grid on;
legend('show')
xticks([0:pi/2:2*pi]);
xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
set(gca,'TickLabelInterpreter','latex');
xlim([0, 2*pi]);

dxt = 2 * cos(t);
dxn = -sin(t);
dwp = dxt + mu * kn / kt * dxn;
dwm = dxt - mu * kn / kt * dxn;
figure(400)
plot(t, dwp, 'ko-', LineWidth=2, DisplayName='$\dot w^+$'), hold on, grid on;
plot(t, dwm, 'r*-', LineWidth=2, DisplayName='$\dot w^-$'), hold on, grid on;
legend('show')
xticks([0:pi/2:2*pi]);
xticklabels({'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'})
set(gca,'TickLabelInterpreter','latex');
xlim([0, 2*pi]);

function [fn, fnp, ft, ft_pre, w] = TangentialForces_withPre(xn, mu, xt, kt, N, w0, kn)
    fn = max(kn * xn, 0);
    fnp = kn * xn;
    ft = zeros(N, 1);
    ft_pre = zeros(N, 1);
    w = zeros(N, 1);
    wp = w0;
    for i = 1:N
        ft_pre(i) = kt * (xt(i) - wp);
        if abs(ft_pre(i)) > mu * fn(i)
            ft(i) = sign(ft_pre(i)) * mu * fn(i);
            wp = xt(i) - sign(ft_pre(i)) * mu * fn(i) / kt;
        else
            ft(i) = ft_pre(i);
        end
        w(i) = wp;
    end

end