function J = HBMJACOB(x, Omega, pfunc)
    Na = pfunc.HBM.Na;
    H = pfunc.HBM.H;
    % linear part
    JL = Jlinear(Omega, pfunc);

    % non-linear part
    JNL = zeros(size(JL));
    xc = x((2 * H + 1) * Na + 1:end);

    % numerical jacobien
    % dGdx = finite_diff_jac(@(x) fftgx(x, pfunc), xc);
    
    % analytical jacobien
    [dGdx, JNLt, Ft, wt, Mft] = JNL_Analytical(xc, pfunc);

    JNL((2 * H + 1) * Na + 1:end, (2 * H + 1) * Na + 1:end) = dGdx;

    J = JNL + JL;
end

