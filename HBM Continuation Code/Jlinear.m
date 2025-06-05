function JLinear = Jlinear(Omega,pfunc)
    
    H = pfunc.HBM.H;
    Dx = 3 * pfunc.HBM.Nx; % dof of contact part
    Na = pfunc.HBM.Na;
    xi = pfunc.HBM.xi;
    Kxx = pfunc.CB_MK.Kxx;
    Kaa = pfunc.CB_MK.Kaa; % conlumn vector
    Kaa = diag(Kaa); % diagnal matrix
    Mxx = pfunc.CB_MK.Mxx;
    Max = pfunc.CB_MK.Max;
    Idx = pfunc.HBM.Idx;
    TAA = zeros((2 * H + 1) * Na);
    TAA(1:Na,1:Na) = Kaa;
    for k = 1:H
        TAA1 = Kaa - k^2 * Omega^2 * eye(size(Kaa));
        TAA2 = k * Omega * xi * Kaa;
        TAA((2 * k - 1) * Na + 1 : (2 * k + 1) * Na, (2 * k - 1) * Na + 1 : (2 * k + 1) * Na) ...
            = [TAA1, TAA2; -TAA2', TAA1];
    end

    TAX = zeros((2 * H + 1) * Na, (2 * H + 1) * Dx);
    for k = 1:H
        TAX1 = - k^2 * Omega^2 * Max;
        TAX2 = zeros(size(TAX1));
        TAX((2 * k - 1) * Na + 1 : (2 * k + 1) * Na, (2 * k - 1) * Dx + 1 : (2 * k + 1) * Dx) ...
            = [TAX1, TAX2; TAX2, TAX1];
    end

    

    TXX = zeros((2 * H + 1) * Dx);
    TXX(1:Dx,1:Dx) = Kxx; 
    for k = 1:H
        TXX1 = Kxx - k^2 * Omega^2 * Mxx;
        TXX2 = k * Omega * xi * Kxx;
        TXX((2 * k - 1) * Dx + 1 : (2 * k + 1) * Dx, (2 * k - 1) * Dx + 1 : (2 * k + 1) * Dx) ...
            = [TXX1, TXX2; -TXX2, TXX1];
    end

    T = [TAA,TAX; TAX',TXX];
    JLinear = T(Idx, Idx); % in final order
end

