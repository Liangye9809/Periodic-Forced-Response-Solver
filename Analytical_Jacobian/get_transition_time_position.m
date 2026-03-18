function [C_ps, C_ss, C_gs] = get_transition_time_position(a, c)
% a is stick condition each time instant, 1 is stick, 0 is non-stick, column vector
% b is slip condition each time instant, 1 or -1 is slip, 0 is non-slip, column vector
% c is contact condition each time instant, 1 is contact, 0 is gap, column vector
% C_ss stick after slip time range
% C_gs stick after gap time range
% C_sp slip time range in plus sign
% C_sm slip time range in minus sign
    N = size(a, 1);
    C_ps = 0; 
    C_ss = - N .* ones(10, 2); n_ss = 1;
    C_gs = - N .* ones(10, 2); n_gs = 1;
    % find stick time range, consider cycle is [1, N]
    if ismember(0, a)
        keep = 0;
        i_stop = N;
        for i = 1:N
            if i > i_stop
                break;
            end
            if i == 1 && a(1) ~=0 % if begin with stick, means transition time located previously 
                for j = 1:N % backwards loop to find slip to stick transition
                    k = N - j + 1;
                    if a(k) == 0
                        i_stop = k;
                        p_plus = mod(k, N) + 1; % t+ transition position, stick
                        p_minus = mod(k - 1, N) + 1; % t- transition position, non-stick
                        p_mid = 0.5 * (p_plus + p_minus); % middle transition position.
                        if k == N
                            p_mid = N - 0.5;
                        end
                        keep = 1;
                        if c(k) == 0
                            C_gs(n_gs, 1) = p_plus; % consider cycle is [1, N]
                        else
                            C_ss(n_ss, 1) = p_plus; % consider cycle is [1, N]
                        end
                        break;
                    end
                end
            elseif a(i) ~= 0 && keep == 0
                p_plus = i; % t+ transition position.
                p_minus = i - 1; % t- transition position.
                p_mid = 0.5 * (p_plus + p_minus); % middle transition position.
                keep = 1;
                if c(p_minus) == 0
                    C_gs(n_gs, 1) = p_plus;
                else
                    C_ss(n_ss, 1) = p_plus;
                end
            end
        
            if a(i) == 0 && keep == 1 % stick end
                keep = 0;
                p_plus = i; % t+ transition position, non-stick.
                p_minus = i - 1; % t- transition position, stick.
                p_mid = 0.5 * (p_plus + p_minus); % middle transition position.
                if C_gs(n_gs, 1) > 0
                    C_gs(n_gs, 2) = p_minus;
                    n_gs = n_gs + 1;
                else
                    C_ss(n_ss, 2) = p_minus;
                    n_ss = n_ss + 1;
                end
            end
    
        end
    else
        C_ps = 1; % pure stick case 
    end

   

   
    if C_ss(n_ss, 1) > 0 && C_ss(n_ss, 2) < 0
        C_ss(n_ss, 2) = N;
    end

    if C_gs(n_gs, 1) > 0 && C_gs(n_gs, 2) < 0
        C_gs(n_gs, 2) = N;
    end

   

    
end
