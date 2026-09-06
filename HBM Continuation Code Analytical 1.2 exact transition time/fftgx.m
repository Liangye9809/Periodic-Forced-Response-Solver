% x are in frequency domain 
% x = [a¹0,a¹1,b¹1,a¹2,b¹2,...,a¹H,b¹H,  a²0,a²1,b²1,a²2,b²2,...,a²H,b²H,...]'

%% original structure
function [F, w, flag] = fftgx(xct, pfunc) % x(t) = E*X
    
    EH = pfunc.HBM.EH;
     N = pfunc.HBM.N;
    
    gxp = pfunc.static.preload.gxp;
    
       kn = pfunc.fc.kn;
      xn0 = pfunc.fc.xn0;
       mu = pfunc.fc.mu;
       kt = pfunc.fc.kt;
     w_in = pfunc.fc.w;
    nloop = pfunc.fc.nloop;

    % xct = xt + xp';

    [Fti, wi, flag] = g(xct, kn, xn0, mu, kt, w_in, nloop); 
    w = wi(1:2, :, end);
    % Calculate Segments
    segments_All = get_all_segments(flag(:, :, end - N + 1:end), xct, dxct, kt, kn, mu, Fti(end - N + 1:end, :), H);
    
    Ft = Fti(end - N + 1:end, :) - gxp'; % the last periods
    
    % TestF = [TestF, Fnt ./ Fti(end - N + 1:end, 3:3:end)];

    hndn = EH * Ft;
    F = hndn(:);
    
end


