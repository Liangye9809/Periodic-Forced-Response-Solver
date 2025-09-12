% x are in frequency domain 
% x = [a¹0,a¹1,b¹1,a¹2,b¹2,...,a¹H,b¹H,  a²0,a²1,b²1,a²2,b²2,...,a²H,b²H,...]'

%% original structure
function F = fftgx(x, pfunc) % x(t) = E*X
    E = pfunc.HBM.E;
    EH = pfunc.HBM.EH;
    xp = pfunc.static.preload.xp;
    gxp = pfunc.static.preload.gxp;
    n = size(E, 2); % 2H+1
    a = size(x, 1) / n; % number of DOF
    X = zeros(n, a);
    for i = 1:a
        r1 = n * (i - 1) + 1;
        r2 = n * i;
        X(:,i) = x(r1:r2); % reorder in dofs in column
    end
    xt = E * X; 
    % gxstruct = g(xt + xp', pfunc.fc); % call matlab function
    % fc.kn = pfunc.fc.kn;
    % fc.xn0 = pfunc.fc.xn0;
    % fc.mu = pfunc.fc.mu;
    % fc.kt = pfunc.fc.kt;
    % fc.w = pfunc.fc.w;
    gxstruct = g_mex(xt + xp', pfunc.fc); % call mex function
    
    gt = gxstruct.F - gxp';
    hndn = EH * gt;
    % save Ff_omega1.2.mat hndn
    F.F = hndn(:);
    F.w = gxstruct.w;
end

%% same logic with Julia code from Javier
% function F = fftgx(X, pfunc) % x(t) = E*X
%     xp = pfunc.static.preload.xp;
%     gxp = pfunc.static.preload.gxp;
%     H = pfunc.HBM.H;
%     E = pfunc.HBM.E;
%     EH = pfunc.HBM.EH;
%     Nc = pfunc.HBM.Nc;
%     G = zeros(size(X));
%     fc = pfunc.fc;
%     w = zeros(2, Nc);
%     for i = 1:Nc
%         r1 = 1 + 3*(i - 1) * (2*H + 1);
%         r2 = 3*i * (2*H + 1);
%         r3 = 3*(i - 1) + 1;
%         r4 = 3*i;
%         fc_i = struct('w', fc.w(:, i), 'kn', fc.kn(i), 'xn0', fc.xn0(i), ...
%                           'mu', fc.mu(:, i), 'kt', fc.kt(:, i));
%         Gstruct = force_fourier(X(r1:r2), xp(r3:r4), gxp(r3:r4), H, E, EH, ...
%                                 fc_i);
%         G(r1:r2) = Gstruct.F;
%         w(:, i) = Gstruct.w;
%         % if norm(Gstruct.w) > 0
%         %     stop = 1;
%         % end
%     end
% 
%     F.F = G;
%     F.w = w;
% end