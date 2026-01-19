function JNL = HBMJACOB_analytical(dG, H)

JNL = zeros(2 * H + 1);

JNL(:, 1) = 0.5 .* dG(1:(2 * H + 1));
JNL(1, 2:end) = dG(2:(2 * H + 1));

dGp = zeros(size(dG, 1) + 1, 1); % include sin0
dGp(1) = dG(1);
dGp(3:end) = dG(2:end);
for k = 1:H
    for j = 1:H
        
        a = j + k;
        b = abs(k - j);
        JNL(2 * k, 2 * j) = 0.5 * (dGp(2 * a + 1) + dGp(2 * b + 1));
        JNL(2 * k + 1, 2 * j) = 0.5 * (dGp(2 * a + 2) + sign(k - j) * dGp(2 * b + 2));
        JNL(2 * k, 2 * j + 1) = 0.5 * (dGp(2 * a + 2) - sign(k - j) * dGp(2 * b + 2));
        JNL(2 * k + 1, 2 * j + 1) = -0.5 * (dGp(2 * a + 1) - dGp(2 * b + 1));

    end
end

    
end

