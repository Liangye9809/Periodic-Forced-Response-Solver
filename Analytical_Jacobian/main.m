clear
clc
%% for single dof of variable
eps = zeros(11, 10);
h = 10^(-8);
order = 1;
h_con = zeros(11, 10);
for iH = 1:11
    for iN = 1:12
        % N = 2 ^ (iN + 4);
        N = 256;
        H = iH;
        ih = iN + 0;
        h = 10^(-ih);
        dt = 2 * pi / N;
        t = 0:dt:(2 * pi - dt);
        xt = exp(sin(t))';
        dgt = dg(xt);
        % plot(t,x);
        % plot(t,dgt);
        
        % Numerical Jacobian
        [E, EH] = HBM.fft_matrices(N, H);
        X = EH * xt;
        JNL_N = finite_diff_jac(@(x) fftgx(x, E, EH), X, h, order);
        
        % Analytiacl Jacobian
        [E, EH] = HBM.fft_matrices(N, 2 * H);
        dGmn = EH * dgt;
        JNL_A = HBMJACOB_analytical(dGmn, H);

        eps(iH, iN) = norm(JNL_A - JNL_N) / norm(JNL_A);
        h_con(iH, iN) = h;
    end
end

%% for multi dof of variables
clear
clc
eps = zeros(11, 10);
h = 10^(-8);
order = 2;
h_con = zeros(11, 10);
N = 256;
dt = 2 * pi / N;
t = 0:dt:(2 * pi - dt);
xt = fxt(t);
dgt = dg(xt);
% plot(t,x);
% plot(t,dgt);
for iH = 1:11
    for iN = 1:12
        
        H = iH;
        ih = iN + 0;
        h = 10^(-ih);

        % Numerical Jacobian
        [E, EH] = HBM.fft_matrices(N, H);
        X = EH * xt;
        X = X(:);
        % [F, xtF] = fftgx(X, E, EH); % for testing
        JNL_N = finite_diff_jac(@(x) fftgx(x, E, EH), X, h, order);
        
        % Analytiacl Jacobian
        [E, EH] = HBM.fft_matrices(N, 2 * H);
        dGmn = []; % initialize the size
        for iG = 1:3
            for jG = 1:3
                dGtemp = [];
                dGtemp(:, 1) = dgt(iG, jG, :);
                dGmn(iG, jG, :) = EH * dGtemp;
            end
        end
        JNL_A = HBMJACOB_analytical_muilti(dGmn, H);

        eps(iH, iN) = norm(JNL_A - JNL_N) / norm(JNL_A);
        h_con(iH, iN) = h;
    end
end