function Fnt_out = ScaleFn(Fnt_in, xn)
    Fnt_out = Fnt_in;
    [N, Nx] = size(xn);
    s = sign(xn);
    for i = 1:Nx
        diffs = [diff(s(:, i)); s(1, i) - s(end, i)];

        trans_contact = find(diffs == -2); % contact to gap
        for j = 1:size(trans_contact, 1)
            i_m = trans_contact(j);
            i_p = mod(i_m, N) + 1; % next point
            c = 0.5 * (1 + xn(i_m, i) / (xn(i_m, i) - xn(i_p, i)));
            Fnt_out(i_m, i) = Fnt_in(i_m, i) * c;
        end

        trans_gap = find(diffs == 2); % gap to contact
        for j = 1:size(trans_gap, 1)
            i_m = trans_gap(j);
            i_p = mod(i_m, N) + 1; % next point
            c = 0.5 * (2 - xn(i_m, i) / (xn(i_m, i) - xn(i_p, i)));
            Fnt_out(i_p, i) = Fnt_in(i_p, i) * c;
        end

        
    end

end