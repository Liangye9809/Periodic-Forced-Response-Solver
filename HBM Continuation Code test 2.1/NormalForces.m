function FN = NormalForces(xn, kn, xn0)
    u = xn - xn0;
    FN = max(0, kn .* u);
end