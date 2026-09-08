function J = HBMJACOB(pfunc, JL, segments_all)
    Na = pfunc.HBM.Na;
    H = pfunc.HBM.H;
    kt = pfunc.fc.kt;
    kn = pfunc.fc.kn;
    mu = pfunc.fc.mu;
    % non-linear part
    JNL = zeros(size(JL));
    
    dGdx = JNL_Analytical(segments_all, H, kt, kn, mu);

    % dGdx = JNL_Analytical_pre(xct, flag(:, :, end - N + 1:end), H, N, kt, kn, mu);

    JNL((2 * H + 1) * Na + 1:end, (2 * H + 1) * Na + 1:end) = dGdx;

    J = JNL + JL;
end

