function Ft_out = ScaleFt(Ft_in, xt) % all the dofs of xt, and Ft
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
                
            Ft_out(i_p, 3 * i - 2) = Ft_in(i_p, 3 * i - 2) * cn; % Ft1
            Ft_out(i_p, 3 * i - 1) = Ft_in(i_p, 3 * i - 1) * cn; % Ft2

        end

        
    end

end