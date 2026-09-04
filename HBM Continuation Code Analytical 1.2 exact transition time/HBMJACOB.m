function J = HBMJACOB(xct, pfunc, JL, flag)
    Na = pfunc.HBM.Na;
    H = pfunc.HBM.H;
    N = pfunc.HBM.N;
    kt = pfunc.fc.kt;
    kn = pfunc.fc.kn;
    mu = pfunc.fc.mu;
    % non-linear part
    JNL = zeros(size(JL));
    
    dGdx = JNL_Analytical(xct, flag(:, :, end - N + 1:end), H, N, kt, kn, mu);

    JNL((2 * H + 1) * Na + 1:end, (2 * H + 1) * Na + 1:end) = dGdx;

    J = JNL + JL;
end

