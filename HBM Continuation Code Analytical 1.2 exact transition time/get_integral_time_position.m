function segments = get_integral_time_position(flag, xt, xn, dxt, dxn, kt, kn, mu, ft, H)
    % FIND_CYCLIC_SEGMENTS Find cyclic time ranges for each segment in flag vector
    % Time axis: t(i) = (i-1)/N * 2pi, i = 1..N
    % Boundaries are at midpoints between transitions
    
    flag = flag(:);
    N = length(flag);
    
    % --- Find transition points (cyclic) ---
    % diff detects where value changes; also check wrap-around
    diffs = [diff(flag); flag(1) - flag(end)];  % length N, last element is wrap
    trans_idx = find(diffs ~= 0);               % indices WHERE change happens
    % Boundary midpoint index (fractional): between trans_idx and trans_idx+1 (cyclic)
    n_trans = length(trans_idx);
    
    if n_trans == 0
        % fprintf('Signal is constant: value = %d over [0, 2pi]\n', flag(1));
        segments = struct('value', flag(1), 't_start', 0, 't_end', 2*pi, 'w', xt(end) - ft(end) / kt);
        return;
    end
    
   
    
    % --- Build segments ---
    % Each segment starts at one boundary and ends at the next
    % Value of segment = flag at the first index AFTER the boundary
    % dxdn = dxt / dxn in gap to stick transition time instant, for Jacobian Calculation
    % Mw and vw = integral of Fourier function, for Jacobian and function
    
    segments = struct('value', {}, 't_start', {}, 't_end', {}, 'dxdn', {}, 'MW', {}, 'vw', {}, 'w', {});
    i = 1;
    for k = 1:n_trans
        n_m = trans_idx(k); % transition point left
        n_p = mod(n_m, N) + 1; % transition point right
        Tstate = 10 * flag(n_m) + flag(n_p);
        switch Tstate
            case {1, -1, 2, 10, -10} % gap to contact, slip to gap
                fdt = xn(n_m) / (xn(n_m) - xn(n_p));
                tau_n = n_m + fdt;
                tau = (tau_n - 1) / N * 2 * pi;
                segments(i).value = flag(n_p);
                segments(i).t_start = tau;
                if Tstate == 2 % gap to stick
                    dxt_tau = dxt(n_m) + fdt * (dxt(n_p) - dxt(n_m));
                    dxn_tau = dxn(n_m) + fdt * (dxn(n_p) - dxn(n_m));
                    segments(i).dxdn = dxt_tau / dxn_tau; % for Jacobian
                    segments(i).w = xt(n_p) - ft(n_p) / kt;
                end
                i = i + 1;
            case 20 % stick (to slip) to gap
                f1_m = ft(n_m);
                f1_p = f1_m + kt * (xt(n_p) - xt(n_m));
                fdt_gap = xn(n_m) / (xn(n_m) - xn(n_p));
                f_gap = f1_m + (f1_p - f1_m) * fdt_gap;
                tau2_n = n_m + f_gap;
                tau2 = (tau2_n - 1) / N * 2 * pi; % contact-gap
                if f_gap ~= 0 % nonzero force at lift-off implies hidden slip
                    f2_m = sign(f_gap) * mu * kn * xn(n_m);
                    fdt = fdt_gap * (f2_m - f1_m) / (f_gap - f1_m + f2_m);
                    tau1_n = n_m + fdt;
                    tau1 = (tau1_n - 1) / N * 2 * pi; % stick-slip
                    segments(i).value = sign(f_gap);
                    segments(i).t_start = tau1;
                    i = i + 1;
                end
                segments(i).value = flag(n_p);
                segments(i).t_start = tau2;
                i = i + 1;
            case {21, 19} % stick to slip
                f1_m = ft(n_m);
                f2_p = ft(n_p);
                f1_p = f1_m + kt * (xt(n_p) - xt(n_m));
                f2_m = flag(n_p) * mu * kn * xn(n_m);
                fdt = (f2_m - f1_m) / (f2_m - f1_m + f1_p - f2_p);
                tau_n = n_m + fdt;
                tau = (tau_n - 1) / N * 2 * pi;
                segments(i).value = flag(n_p);
                segments(i).t_start = tau;
                i = i + 1;
            case {-8, 12} % slip to stick
                dw_m = dxt(n_m) - flag(n_m) * mu * kn / kt * dxn(n_m);
                dw_p = dxt(n_p) - flag(n_m) * mu * kn / kt * dxn(n_p);
                fdt = dw_m / (dw_m - dw_p);
                tau_n = n_m + fdt;
                tau = (tau_n - 1) / N * 2 * pi;
                segments(i).value = flag(n_p);
                segments(i).t_start = tau;
                segments(i).w = xt(n_p) - ft(n_p) / kt;
                i = i + 1;
            otherwise
                error('exist other transition type not define!')
        end
        
    end

    for k = 1:i-1
        segments(k).t_end = segments(mod(k, i - 1) + 1).t_start;
        if segments(k).t_start > segments(k).t_end
            segments(k).t_end  = segments(k).t_end + 2 * pi;
        end
        [segments(k).MW, segments(k).vw] = fW(segments(k).t_start, segments(k).t_end, H);
        
    end
end
