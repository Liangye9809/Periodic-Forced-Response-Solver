function J = HBMJACOB(pfunc, JL, segments_all, xct, flag, x)
    Na = pfunc.HBM.Na;
    H = pfunc.HBM.H;
    N = pfunc.HBM.N;
    kt = pfunc.fc.kt;
    kn = pfunc.fc.kn;
    mu = pfunc.fc.mu;
    % non-linear part
    JNL = zeros(size(JL));
    
    dGdx = JNL_Analytical(segments_all, H, kt, kn, mu);

    % dGdx = JNL_Analytical_pre(xct, flag(:, :, end - N + 1:end), H, N, kt, kn, mu);

    % numerical
    % xc = x(Na * (2 * H + 1) + 1:end);
    % dGdx = finite_diff_jac(@(x) fftgx(x, xct, pfunc), xc);

    JNL((2 * H + 1) * Na + 1:end, (2 * H + 1) * Na + 1:end) = dGdx;

    J = JNL + JL;
end

