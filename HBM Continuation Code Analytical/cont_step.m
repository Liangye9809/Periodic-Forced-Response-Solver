function [x, omega, tx, tomega, w, k, FlagState] = cont_step(func, jacob, deromega, params)
     x0 = params.cont.x0;
omega_0 = params.cont.omega_0;
     ds = params.cont.ds;
    tx0 = params.cont.tx0;
tomega0 = params.cont.tomega0;
     xp = params.func.static.preload.xp;
   epsx = params.Newton.epsx;
   epsf = params.Newton.epsf;
maxiter = params.Newton.maxiter;
      x = x0 + ds*tx0;
  omega = omega_0 + ds*tomega0;
      z = [x; omega];

Nx = params.func.HBM.Nx;
H  = params.func.HBM.H;
Na = params.func.HBM.Na;
E  = params.func.HBM.E;


FlagState = zeros(Nx, 4);

xct = Fourier_to_Time(x(Na * (2 * H + 1) + 1:end), H, Nx, E);
%% calculate (min(w+) - max(w-)) / 2
pfunc = params.func;
w = get_w_middle(xct + xp', pfunc);
params.func.fc.w = w;
%%
[F, w, JL, flag] = func(x, xct + xp', omega, params.func);

params.func.fc.w = w; % update w
G = [F 
     tx0' * (x - x0) + tomega0 * (omega - omega_0) - ds];
for k = 1:maxiter
    JG = [jacob(xct + xp', params.func, JL, flag) deromega(x, omega, params.func)
          tx0' tomega0];
    
    dz = -JG\G;
    z = z + dz;
    x = z(1:end-1);
    omega = z(end);
    
    xct = Fourier_to_Time(x(Na * (2 * H + 1) + 1:end), H, Nx, E);
    [F, w, JL, flag] = func(x, xct + xp', omega, params.func);
    
    params.func.fc.w = w; % update w
   

    G = [F 
         tx0' * (x - x0) + tomega0 * (omega - omega_0) - ds];
    errorx = norm(dz) / norm(z);
    errorf = norm(G);
    % format short g
    % disp([k z(1) z(2) lambda errorx  errorf])
    if (errorx < epsx) && (errorf < epsf)
        t0 = [tx0; tomega0];

        A = [jacob(xct + xp', params.func, JL, flag) deromega(x, omega, params.func)];
        [Q, R] = qr(A');
        t = Q(:, end);
        t = t * sign(t0'*t);
        tx = t(1:end-1);
        tomega = t(end);
        format short g
        disp([params.cont.step  z(1) z(2) omega  errorx  errorf   k])
       
        %% output slip state
        FlagState = get_slip_state(flag); % [2, 1, -1, 0] = [stick, slipPlus, slipMinus, gap] [Nx * 4]
        
        %%
        return
    end
end


end

function xct = Fourier_to_Time(Xc, H, Nx, E)

    for i = 1:3 * Nx
        r1 = (2 * H + 1) * (i - 1) + 1;
        r2 = (2 * H + 1) * i;
        X(:,i) = Xc(r1:r2); % reorder in dofs in column
    end
    xct = E * X; 

end

function w = get_w_middle(xct, pfunc)
    N = pfunc.HBM.N;
    Nx = pfunc.HBM.Nx;
    kn = pfunc.fc.kn;
    kt = pfunc.fc.kt;
    w_in = pfunc.fc.w;
    mu = pfunc.fc.mu;

    w = zeros(2, Nx);
    xtt = zeros(N, 2);
    xtn = zeros(N, 1);
    for i = 1:Nx
        xtt(:, 1) = xct(:, 3 * i - 2);
        xtt(:, 2) = xct(:, 3 * i - 1);
        xtn = xct(:, 3 * i);
        for j = 1:2
            w_plus = xtt(:, j) + mu(j, i) * kn(i) / kt(j, i) * max(xtn, 0);
            w_minus = xtt(:, j) - mu(j, i) * kn(i) / kt(j, i) * max(xtn, 0);
            if min(w_plus) >= max(w_minus)
                w(j, i) = 0.5 * (min(w_plus) + max(w_minus));
            else
                w(j, i) = w_in(j, i);
            end
            % t = [0:N-1]';
            % figure;
            % plot(t, w_plus, 'k-'), hold on;
            % plot(t, w_minus, 'b--'), grid on;
            % plot(t, w(j, i) * ones(N, 1), 'r-');
        end
    end

        
end

function FlagState = get_slip_state(flag)
    Nx = size(flag, 2);
    FlagState = zeros(Nx, 4);
    for i = 1:Nx
        FlagState(i, :) = ismember([2, 1, -1, 0], flag(:, i, :));
    end

end









