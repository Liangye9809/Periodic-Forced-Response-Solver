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
        N = 32;
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
        dG = EH * dgt;
        JNL_A = HBMJACOB_analytical(dG, H);

        eps(iH, iN) = norm(JNL_A - JNL_N) / norm(JNL_A);
        h_con(iH, iN) = h;
    end
end