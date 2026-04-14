function J = HBMJACOB(x, pfunc, JL, flag)
    Na = pfunc.HBM.Na;
    H = pfunc.HBM.H;
    N = pfunc.HBM.N;
    % non-linear part
    JNL = zeros(size(JL));
    
    dGdx = JNL_Analytical(flag(:, :, end - N + 1:end));

    JNL((2 * H + 1) * Na + 1:end, (2 * H + 1) * Na + 1:end) = dGdx;

    J = JNL + JL;
end

