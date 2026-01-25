function JNL = HBMJACOB_analytical_muilti(dGmn, H)

    [m, n, l] = size(dGmn);
    if m ~= n
        error('the size of column and row is not equal');
    end
    if (l - 1) / 4 ~= H
        error('the size of Harmonics is not equal');
    end
    
    Nx = m;
    JNL = size(Nx * (2 * H + 1), Nx * (2 * H + 1));

    
    for i = 1:Nx
        for j = 1:Nx
            inda1 = (i - 1) * (2 * H + 1) + 1;
            inda2 = i * (2 * H + 1);
            indb1 = (j - 1) * (2 * H + 1) + 1;
            indb2 = j * (2 * H + 1);
            dGtemp(:, 1) = dGmn(i, j, :);
            JNL(inda1:inda2, indb1:indb2) = HBMJACOB_analytical(dGtemp, H);
        end
    end
end

  