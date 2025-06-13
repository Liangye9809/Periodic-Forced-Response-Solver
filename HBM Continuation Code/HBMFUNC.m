function FUNC = HBMFUNC(x, Omega, pfunc) % x = [a¹0,a¹1,b¹1,a¹2,b¹2,...,a¹H,b¹H,  a²0,a²1,b²1,a²2,b²2,...,a²H,b²H,...]'
    fftfa = pfunc.HBM.fftfa;
    fftfx = pfunc.HBM.fftfx;
    F = [fftfa; fftfx];
    Na = pfunc.HBM.Na; % DOF of elastic part
    H = pfunc.HBM.H;
    JL = Jlinear(Omega, pfunc); % get from outside
    
    xc = x((2 * H + 1) * Na + 1:end);
    G = zeros(size(x));
    Gistruct = fftgx(xc, pfunc);
    Gi = Gistruct.F;
    G((2 * H + 1) * Na + 1:end) = Gi;

    FUNC.F = JL * x + G - F;
    FUNC.w = Gistruct.w;

end

