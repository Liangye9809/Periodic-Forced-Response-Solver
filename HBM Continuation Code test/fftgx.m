% x are in frequency domain 
% x = [a¹0,a¹1,b¹1,a¹2,b¹2,...,a¹H,b¹H,  a²0,a²1,b²1,a²2,b²2,...,a²H,b²H,...]'

%% original structure
function [F, w, flag] = fftgx(x, xct, pfunc) % 


       EH = pfunc.HBM.EH;
        N = pfunc.HBM.N;
      gxp = pfunc.static.preload.gxp;
       xp = pfunc.static.preload.xp;
       kn = pfunc.fc.kn;
      xn0 = pfunc.fc.xn0;
       mu = pfunc.fc.mu;
       kt = pfunc.fc.kt;
     w_in = pfunc.fc.w;
    nloop = pfunc.fc.nloop;

%% inputs are xct, pfunc
    % [Fti, wi, flag] = g(xct + xp', kn, xn0, mu, kt, w_in, nloop); 
    % w = wi(1:2, :, end);
    % 
    % Ft = Fti(N + 1:end, :) - gxp';
    % hndn = EH * Ft;
    % F = hndn(:);
%% inputs are x, pfunc

    E = pfunc.HBM.E;
    n = size(E, 2); % 2H+1
    a = size(x, 1) / n; % number of DOF
    X = zeros(n, a);
    for i = 1:a
        r1 = n * (i - 1) + 1;
        r2 = n * i;
        X(:,i) = x(r1:r2); % reorder in dofs in column
    end
    xt = E * X; 


    [Fti, wi, flag] = g(xt + xp', kn, xn0, mu, kt, w_in, nloop); 
    w = wi(1:2, :, end);

    flag_ = flag(:, :, end - N + 1:end);
    M_fstar = zeros(n, a);
    if sum(ismember([-1, 1, 0], flag_)) > 0
        xct_ = xt + xp';
        % Fnt = ScaleFn(Fti(end - N + 1:end, 3:3:end), xct_(:, 3:3:end)); % pass only normal displacements and normal forces
        % Fti(end - N + 1:end, 3:3:end) = Fnt;
        
        % first slip and stick transition
        [Fti(end - N + 1:end, :), M_fstar] = FFtFactor(Fti(end - N + 1:end, :), xct_, flag_, kt, kn, mu);
    
        % second gap and contact transition
        Fti(end - N + 1:end, :) = ScaleFt(Fti(end - N + 1:end, :), xct_); % pass all the displacements and forces
    end
    Ft = Fti(end - N + 1:end, :) - gxp';
    hndn = EH * Ft + M_fstar;
    F = hndn(:);

    
end


