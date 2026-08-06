function J = HBMJACOB(x, xct, pfunc, JL, flag)
    Na = pfunc.HBM.Na;
    H = pfunc.HBM.H;
    N = pfunc.HBM.N;
    kt = pfunc.fc.kt;
    kn = pfunc.fc.kn;
    mu = pfunc.fc.mu;
    
    % non-linear part
    JNL = zeros(size(JL));
    
    % analytical
    % dGdx = JNL_Analytical(xct, flag(:, :, end - N + 1:end), H, N, kt, kn, mu);

    % numerical jacobien
    xc = x(Na * (2 * H + 1) + 1:end);
    dGdx = finite_diff_jac(@(x) fftgx(x, xct, pfunc), xc);


    JNL((2 * H + 1) * Na + 1:end, (2 * H + 1) * Na + 1:end) = dGdx;

    J = JNL + JL;
end

