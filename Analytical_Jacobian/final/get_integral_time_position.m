function [C_ss, C_gs, C_sp, C_sm] = get_integral_time_position(a, b, c)
% a is stick condition each time instant, 1 is stick, 0 is non-stick, column vector
% b is slip condition each time instant, 1 or -1 is slip, 0 is non-slip, column vector
% c is contact condition each time instant, 1 is contact, 0 is gap, column vector
% C_ss stick after slip time range
% C_gs stick after gap time range
% C_sp slip time range in plus sign
% C_sm slip time range in minus sign
    C_ss = - ones(10, 2); n_ss = 1;
    C_gs = - ones(10, 2); n_gs = 1;
    C_sp = - ones(10, 2); n_sp = 1;
    C_sm = - ones(10, 2); n_sm = 1;
    keep = 0;
    N = size(a, 1);
    % find stick time range, consider cycle is [0, N-1]
    for i = 1:N
        if i == 1 && a(1) ~=0 % if begin with stick, means transition time located previously 
            for j = 1:N % backwards loop to find slip to stick transition
                k = N - j + 1;
                if a(k) == 0
                    p_plus = mod(k, N); % t+ transition instant, stick
                    p_minus = mod(k - 1, N); % t- transition instant, non-stick
                    p_mid = 0.5 * (p_plus + p_minus); % middle transition instant.
                    if k == N
                        p_mid = N - 0.5;
                    end
                    keep = 1;
                    if c(k) == 0
                        C_gs(n_gs, 1) = p_mid; % consider cycle is [0, N-1]
                    else
                        C_ss(n_ss, 1) = p_mid; % consider cycle is [0, N-1]
                    end
                    break;
                end
            end
        elseif a(i) ~= 0 && keep == 0
            p_plus = i - 1; % t+ transition instant.
            p_minus = i - 2; % t- transition instant.
            p_mid = 0.5 * (p_plus + p_minus); % middle transition instant.
            keep = 1;
            if c(p_minus) == 0
                C_gs(n_gs, 1) = p_mid;
            else
                C_ss(n_ss, 1) = p_mid;
            end
        end
    
        if a(i) == 0 && keep == 1 % stick end
            keep = 0;
            p_plus = i - 1; % t+ transition instant, non-stick.
            p_minus = i - 2; % t- transition instant, stick.
            p_mid = 0.5 * (p_plus + p_minus); % middle transition instant.
            if C_gs(n_gs, 1) > 0
                C_gs(n_gs, 2) = p_mid;
                n_gs = n_gs + 1;
            else
                C_ss(n_ss, 2) = p_mid;
                n_ss = n_ss + 1;
            end
        end

    end

    % find slip time range
    keep = 0;
    for i = 1:N
        if i == 1 && b(1) ~=0 % if begin with slip, means transition time located previously 
            for j = 1:N % backwards loop to find non-slip to slip transition
                k = N - j + 1;
                if b(k) == 0
                    p_plus = mod(k, N); % t+ transition instant, slip. 
                    p_minus = mod(k - 1, N); % t- transition instant, non-slip
                    p_mid = 0.5 * (p_plus + p_minus); % middle transition instant.
                    if k == N
                        p_mid = N - 0.5;
                    end
                    keep = 1;
                    if b(k) == 1
                        C_sp(n_sp, 1) = p_mid;
                    else
                        C_sm(n_sm, 1) = p_mid;
                    end
                    break;
                end
            end
        elseif b(i) ~= 0 && keep == 0
            p_plus = i - 1; % t+ transition instant, slip.
            p_minus = i - 2; % t- transition instant, non-slip.
            p_mid = 0.5 * (p_plus + p_minus); % middle transition instant.
            keep = 1;
            if b(i) == 1
                C_sp(n_sp, 1) = p_mid;
            else
                C_sm(n_sm, 1) = p_mid;
            end
        end
    
        if b(i) == 0 && keep == 1 % slip end
            keep = 0;
            p_plus = i - 1; % t+ transition instant, non-slip.
            p_minus = i - 2; % t- transition instant, slip.
            p_mid = 0.5 * (p_plus + p_minus); % middle transition instant.
            if C_sp(n_sp, 1) > 0
                C_sp(n_sp, 2) = p_mid;
                n_sp = n_sp + 1;
            else
                C_sm(n_sm, 2) = p_mid;
                n_sm = n_sm + 1;
            end
        end

    end

    if C_sp(n_sp, 1) > 0 && C_sp(n_sp, 2) == -1
        C_sp(n_sp, 2) = N - 0.5;
    end

    if C_ss(n_ss, 1) > 0 && C_ss(n_ss, 2) == -1
        C_ss(n_ss, 2) = N - 0.5;
    end

    if C_gs(n_gs, 1) > 0 && C_gs(n_gs, 2) == -1
        C_gs(n_gs, 2) = N - 0.5;
    end

    if C_sm(n_sm, 1) > 0 && C_sm(n_sm, 2) == -1
        C_sm(n_sm, 2) = N - 0.5;
    end

    % C_ss = C_ss ./ N .* (2 * pi);
    % C_gs = C_gs ./ N .* (2 * pi);
    % C_sp = C_sp ./ N .* (2 * pi);
    % C_sm = C_sm ./ N .* (2 * pi);
end
