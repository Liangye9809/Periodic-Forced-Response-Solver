function Ft_out = ScaleFt(Ft_in, xt, flag, w) % all the dofs of xt, and Ft
    Ft_out = Ft_in;
    [N, Nx] = size(xt);
    Nx = Nx / 3;
    xn = xt(:, 3:3:end);
    s = sign(xn);
    for i = 1:Nx
        diffs = [diff(s(:, i)); s(1, i) - s(end, i)];

        trans_contact = find(diffs == -2); % contact to gap
        for j = 1:size(trans_contact, 1)
            % for Fn
            i_m = trans_contact(j);
            i_p = mod(i_m, N) + 1; % next point
            cn = 0.5 * (1 + xn(i_m, i) / (xn(i_m, i) - xn(i_p, i)));
            Ft_out(i_m, 3 * i) = Ft_in(i_m, 3 * i) * cn;

            % For Ft
            % if flag(1, i, i_m) == 2 % Ft1 is stick
            %     ct = 0.5 * (1 + (xt(i_m, 3 * i - 2) - w(1, i, i_m)) / (1e-16 + xt(i_m, 3 * i - 2) - xt(i_p, 3 * i - 2)));
            % else % slip
            %     ct = cn;
            % end
            % Ft_out(i_m, 3 * i - 2) = Ft_in(i_m, 3 * i - 2) * ct; % Ft1
            % 
            % if flag(2, i, i_m) == 2 % Ft2 is stick
            %     ct = 0.5 * (1 + (xt(i_m, 3 * i - 1) -w(2, i, i_m))/ (1e-16 + xt(i_m, 3 * i - 1) - xt(i_p, 3 * i - 1)));
            % else
            %     ct = cn;
            % end
            % Ft_out(i_m, 3 * i - 1) = Ft_in(i_m, 3 * i - 1) * ct; % Ft2

            Ft_out(i_m, 3 * i - 2) = Ft_in(i_m, 3 * i - 2) * cn; % Ft1
            Ft_out(i_m, 3 * i - 1) = Ft_in(i_m, 3 * i - 1) * cn; % Ft2
        end

        trans_gap = find(diffs == 2); % gap to contact
        for j = 1:size(trans_gap, 1)
            % For Fn
            i_m = trans_gap(j);
            i_p = mod(i_m, N) + 1; % next point
            cn = 0.5 * (2 - xn(i_m, i) / (xn(i_m, i) - xn(i_p, i)));
            Ft_out(i_p, 3 * i) = Ft_in(i_p, 3 * i) * cn;

            % For Ft
            % if flag(1, i, i_p) == 2 % Ft1 is stick
            %     ct = 0.5 * (2 - xt(i_m, 3 * i - 2) / (1e-16 + xt(i_m, 3 * i - 2) - xt(i_p, 3 * i - 2)));
            % else
            %     ct = cn; % slip
            % end
            % Ft_out(i_p, 3 * i - 2) = Ft_in(i_p, 3 * i - 2) * ct; % Ft1
            % 
            % if flag(2, i, i_p) == 2 % Ft2 is stick
            %     ct = 0.5 * (2 - xt(i_m, 3 * i - 1) / (1e-16 + xt(i_m, 3 * i - 1) - xt(i_p, 3 * i - 1)));
            % else
            %     ct = cn;
            % end
            % Ft_out(i_p, 3 * i - 1) = Ft_in(i_p, 3 * i - 1) * ct; % Ft2
                
            Ft_out(i_p, 3 * i - 2) = Ft_in(i_p, 3 * i - 2) * cn; % Ft1
            Ft_out(i_p, 3 * i - 1) = Ft_in(i_p, 3 * i - 1) * cn; % Ft2

        end

        
    end

end