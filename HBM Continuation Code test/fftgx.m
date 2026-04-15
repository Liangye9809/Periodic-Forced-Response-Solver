% x are in frequency domain 
% x = [a¹0,a¹1,b¹1,a¹2,b¹2,...,a¹H,b¹H,  a²0,a²1,b²1,a²2,b²2,...,a²H,b²H,...]'

%% original structure
function [F, w, flag] = fftgx(xct, pfunc) % x(t) = E*X
    % E = pfunc.HBM.E;
    EH = pfunc.HBM.EH;
     N = pfunc.HBM.N;
    % xp = pfunc.static.preload.xp;
    gxp = pfunc.static.preload.gxp;
    % n = size(E, 2); % 2H+1
    % a = size(x, 1) / n; % number of DOF
    % X = zeros(n, a);
    % for i = 1:a
    %     r1 = n * (i - 1) + 1;
    %     r2 = n * i;
    %     X(:,i) = x(r1:r2); % reorder in dofs in column
    % end
    % xt = E * X; 
    % 
       kn = pfunc.fc.kn;
      xn0 = pfunc.fc.xn0;
       mu = pfunc.fc.mu;
       kt = pfunc.fc.kt;
     w_in = pfunc.fc.w;
    nloop = pfunc.fc.nloop;

    % xct = xt + xp';

    [Fti, wi, flag] = g(xct, kn, xn0, mu, kt, w_in, nloop); 
    w = wi(1:2, :, end);
    
    Ft = Fti(N + 1:end, :) - gxp';
    hndn = EH * Ft;
    F = hndn(:);
    
end


