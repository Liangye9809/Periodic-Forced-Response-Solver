function Jomega = HBMJOmega(x,Omega,p)
    H = p.HBM.H;
    Dx = 3 * p.HBM.Nx; % dof of contact part
    Na = p.HBM.Na;
    xi = p.HBM.xi;
    Kxx = p.CB_MK.Kxx;
    Kaa = p.CB_MK.Kaa;
    Kaa = diag(Kaa);
    Mxx = p.CB_MK.Mxx;
    Max = p.CB_MK.Max;
    Idx = p.HBM.Idx;
    dTAA = zeros((2 * H + 1) * Na);
    for k = 1:H
        dTAA1 = - 2 * k^2 * Omega * eye(size(Kaa));
        dTAA2 = k * xi * Kaa;
        dTAA((2 * k - 1) * Na + 1 : (2 * k + 1) * Na, (2 * k - 1) * Na + 1 : (2 * k + 1) * Na) ...
            = [dTAA1,dTAA2; -dTAA2', dTAA1];
    end

    dTAX = zeros((2 * H + 1) * Na, (2 * H + 1) * Dx);
    for k = 1:H
        dTAX1 = - 2 * k^2 * Omega * Max;
        dTAX2 = zeros(size(dTAX1));
        dTAX((2 * k - 1) * Na + 1 : (2 * k + 1) * Na, (2 * k - 1) * Dx + 1 : (2 * k + 1) * Dx) ...
            = [dTAX1,dTAX2; dTAX2, dTAX1];
    end

    dTXX = zeros((2 * H + 1) * Dx);
    for k = 1:H
        dTXX1 = - 2 * k^2 * Omega * Mxx;
        dTXX2 = k * xi * Kxx;
        dTXX((2 * k - 1) * Dx + 1 : (2 * k + 1) * Dx, (2 * k - 1) * Dx + 1 : (2 * k + 1) * Dx) ...
            = [dTXX1,dTXX2; -dTXX2, dTXX1];
    end

    dT = [dTAA,dTAX; dTAX',dTXX];
    dT = dT(Idx,Idx); % in final order

    Jomega = dT * x;
end