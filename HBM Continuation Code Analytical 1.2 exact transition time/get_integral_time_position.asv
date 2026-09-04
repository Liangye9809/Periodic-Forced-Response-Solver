function segments = get_integral_time_position(flag, xt, xn, dxt, dxn, kt, kn, mu, ft)
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
        segments = struct('value', flag(1), 't_start', 0, 't_end', 2*pi, 'index_start', 1);
        return;
    end
    
   
    
    % --- Build segments ---
    % Each segment starts at one boundary and ends at the next
    % Value of segment = flag at the first index AFTER the boundary
    
    segments = struct('value', {}, 't_start', {}, 't_end', {}, 'index_start', {});
    
    for k = 1:n_trans
        n_m = trans_idx(k); % transition point left
        n_p = mod(n_m, N) + 1; % transition point right
        Tstate = 10 * flag(n_m) + flag(n_p);
        switch Tstate
            case {1, -1, 2, 10, -10} % gap to contact, slip to contact
                tau_n = n_m + xn(n_m) / (xn(n_m) - xn(n_p));
                tau = (tau_n - 1) / N * 2 * pi;
            case 20 % stick (to slip) to gap
                f1_m = ft(i_m);
                f2_p = ft(i_p);
                f1_p = f1_m + kt * (xt(i_p) - xt(i_m));
                fdt_gap = xn(i_m) / (xn(i_m) - xn(i_p));
                f_gap = f1_m + (f1_p - f1_m) * fdt_gap;
                tau2_n = n_m + f_gap;
                tau2 = (tau2_n - 1) / N * 2 * pi;
                if f_gap ~= 0 % nonzero force at lift-off implies hidden slip
                    f2_m = sign(f_gap) * mu * kn * xn(i_m);
                    fdt = fdt_gap * (f2_m - f1_m) / (f_gap - f1_m + f2_m);
                    tau1_n = n_m + fdt;
                    tau1 = (tau1_n - 1) / N * 2 * pi;
                end
            case {21, 19} % stick to slip
                f1_m = ft(i_m);
                f2_p = ft(i_p);
                f1_p = f1_m + kt * (xt(i_p) - xt(i_m));
                f2_m = flag(n_p) * mu * kn * xn(i_m);
                fdt = (f2_m - f1_m) / (f2_m - f1_m + f1_p - f2_p);
                tau_n = n_m + fdt;
                tau = (tau_n - 1) / N * 2 * pi;
            case {-8, 12} % slip to stick
                dw_m = dxt(n_m) - flag(n_m) * mu * kn / kt * dxn(n_m);
                dw_p = dxt(n_p) - flag(n_m) * mu * kn / kt * dxn(n_p);
                fdt = dw_m / (dw_m - dw_p);
                tau_n = n_m + fdt;
                tau = (tau_n - 1) / N * 2 * pi;
            otherwise
                error('exist other transition type not define!')
        end
        b_start = ;                          % start boundary (index)
        b_end   = ;        % next boundary (index, cyclic)
        
        % First sample index after b_start
        first_idx = mod(trans_idx(k), N) + 1;               % 1-based, cyclic
        val = flag(first_idx);
        
        % Convert boundary indices to time: t = (idx - 1) / N * 2pi
        t_start = (b_start - 1) / N * 2 * pi;
        t_end   = (b_end   - 1) / N * 2 * pi;
        
        % Handle cyclic wrap: if t_end <= t_start, it wraps around
        if t_end <= t_start
            t_end = t_end + 2 * pi;
        end
        
        segments(k).value       = val;
        segments(k).t_start     = t_start;
        segments(k).t_end       = t_end;
        segments(k).index_start = mod(b_start - 0.5, N) + 1;
    end
end
