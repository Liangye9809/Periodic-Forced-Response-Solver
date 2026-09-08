% x are in frequency domain 
% x = [a¹0,a¹1,b¹1,a¹2,b¹2,...,a¹H,b¹H,  a²0,a²1,b²1,a²2,b²2,...,a²H,b²H,...]'

%% original structure
function [F, w, flag, segments_All] = fftgx(X, xct, pfunc) % x(t) = E*X
    
     E = pfunc.HBM.E;
    EH = pfunc.HBM.EH;
     N = pfunc.HBM.N;
     H = pfunc.HBM.H;
    Nx = pfunc.HBM.Nx;
    Na = pfunc.HBM.Na;
    gxp = pfunc.static.preload.gxp;
    xp = pfunc.static.preload.xp;
    
       kn = pfunc.fc.kn;
      xn0 = pfunc.fc.xn0;
       mu = pfunc.fc.mu;
       kt = pfunc.fc.kt;
     w_in = pfunc.fc.w;
    nloop = pfunc.fc.nloop;

    
    Xc = X(Na * (2 * H + 1) + 1:end);
    for ixp = 1:3 * Nx
        Xc((2 * H + 1) * (ixp - 1) + 1) = Xc((2 * H + 1) * (ixp - 1) + 1) + 2 * xp(ixp); % consider preload
    end
    dxct = get_dxt(Xc, E, Nx);
    

    [Fti, wi, flag] = g(xct, kn, xn0, mu, kt, w_in, nloop, dxct); 
    w = wi(1:2, :, end);
    % Calculate Segments
    segments_All = get_all_segments(flag(:, :, end - N + 1:end), xct, dxct, kt, kn, mu, Fti(end - N + 1:end, :), H);
    F_Analytical = get_Analytical_F_Fourier(segments_All, gxp', kt, kn, mu, H, Xc);
    
    %% numerical Ft
    Ft = Fti(end - N + 1:end, :) - gxp'; % the last periods
    hndn = EH * Ft;
    F_N = hndn(:);
    
    %% analytical F
    F = F_Analytical;
end


