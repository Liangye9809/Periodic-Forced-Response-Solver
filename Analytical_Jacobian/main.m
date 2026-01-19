clear
clc
%% for single dof of variable
eps = zeros(7, 4);
h = 1e-6;
order = 2;
for iH = 1:11
    for iN = 1:4
        N = 2 ^ (iN + 4);
        H = iH;
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
    end
end