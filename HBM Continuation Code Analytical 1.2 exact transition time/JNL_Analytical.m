function JNL = JNL_Analytical(x, flag, H, N, kt, kn, mu) % x is the size of N*3Nx

    Nx = size(flag, 2);
    JNL = zeros(3 * Nx * (2 * H + 1), 3 * Nx * (2 * H + 1));

    for i = 1:Nx
        flagT1(:, 1) = flag(1, i, :);
        flagT2(:, 1) = flag(2, i, :);
        JNLi = JNL_one_Nx(x(:, 3 * i - 2:3 * i), flagT1, flagT2, H, N, kt(:, i), kn(i), mu(:, i));
        indx1 = 3 * (i - 1) * (2 * H + 1) + 1;
        indx2 = 3 * i * (2 * H + 1);
        JNL(indx1:indx2, indx1:indx2) = JNLi;
    end
            

end


function JNLi = JNL_one_Nx(x, flagT1, flagT2, H, N, kt, kn, mu)

    JNLi = zeros(3 * (2 * H + 1), 3 * (2 * H + 1));

    segmentsT1 = get_integral_time_position(flagT1);
    [dF1dX1, dF1dXn, dFndXn] = get_dFdX(segmentsT1, H, N, kt(1), kn, mu(1), x(:, 1), x(:, 3));
    
    segmentsT2 = get_integral_time_position(flagT2);
    [dF2dX2, dF2dXn, ~] = get_dFdX(segmentsT2, H, N, kt(2), kn, mu(2), x(:, 2), x(:, 3));

    JNLi(1:2 * H + 1, 1:2 * H + 1) = dF1dX1;
    JNLi(1:2 * H + 1, 2 * (2 * H + 1) + 1:end) = dF1dXn;

    JNLi((2 * H + 1) + 1:2 * (2 * H + 1), (2 * H + 1) + 1:2 * (2 * H + 1)) = dF2dX2;
    JNLi((2 * H + 1) + 1:2 * (2 * H + 1), 2 * (2 * H + 1) + 1:end) = dF2dXn;

    JNLi(2 * (2 * H + 1) + 1:end, 2 * (2 * H + 1) + 1:end) = dFndXn;
end


function segments = get_integral_time_position(flag)
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
    
    % Boundary midpoints (fractional 1-based index, cyclic)
    boundary_idx = trans_idx + 0.5;  % midpoint between last-of-old and first-of-new
    
    % --- Build segments ---
    % Each segment starts at one boundary and ends at the next
    % Value of segment = flag at the first index AFTER the boundary
    
    segments = struct('value', {}, 't_start', {}, 't_end', {}, 'index_start', {});
    
    for k = 1:n_trans
        b_start = boundary_idx(k);                          % start boundary (index)
        b_end   = boundary_idx(mod(k, n_trans) + 1);        % next boundary (index, cyclic)
        
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

function [dFdX, dFdXn, dFndXn] = get_dFdX(segmentsT, H, N, kt, kn, mu, xt, xn)
    dFdX   = zeros(2 * H + 1, 2 * H + 1);
    dFdXn  = zeros(2 * H + 1, 2 * H + 1);
    dFndXn = zeros(2 * H + 1, 2 * H + 1);

    k = length(segmentsT);

    if k == 1 % pure stick or whole gap
        if segmentsT(1).value == 2 % pure stick
            dFdX   = kt .* eye(2 * H + 1);
            dFndXn = kn .* eye(2 * H + 1);
            return;
        elseif segmentsT(1).value == 0 % whole gap
            return;
        end
    end

    for i = 1:k
        t1 = segmentsT(i).t_start;
        t2 = segmentsT(i).t_end;
        switch segmentsT(i).value
            case 2 % stick after slip or gap
                [MW, Mw] = fW(t1, t2, H);
                c_vec = c_vector(t1, H);
                dFdX  = dFdX + kt .* MW - Mw * (kt .* c_vec);
                dFndXn = dFndXn + kn .* MW; % contact
                switch segmentsT(mod(i - 2, k) + 1).value

                    case 0 % gap to stick
                        i_p   = segmentsT(i).index_start; % stick start index
                        i_m   = mod(i_p - 2, N) + 1;  % previous index
                        dxdn  = (xt(i_p) - xt(i_m)) / (xn(i_p) - xn(i_m));
                        dFdXn = dFdXn + Mw * (kt .* dxdn .* c_vec);

                    case {-1, 1} % slip to stick
                        dFdXn = dFdXn + Mw * (segmentsT(mod(i - 2, k) + 1).value .* mu .* kn .* c_vec);
                end 

            case {-1, 1} % slip
                [MW, ~] = fW(t1, t2, H);
                dFdXn   = dFdXn + segmentsT(i).value .* mu .* kn .* MW;
                dFndXn  = dFndXn + kn .* MW; % contact
        end

    end

end