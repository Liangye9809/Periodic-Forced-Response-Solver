function F = gg(x1, x2, x3, fc)
    
    
    kn = fc.kn; % stiffness in normal direction
    xn0 = fc.xn0; % initial gap between this contact point
    mu = fc.mu; % friciton coefficient of coulomb friction
    kt = fc.kt; % stiffness of 2 tangential direction in column vector
    w = fc.w;

    N = NormalForces(x3, kn, xn0);
    [T1, w1] = TangentialForces(x1, w(1), kt(1), mu(1), N);
    [T2, w2] = TangentialForces(x2, w(2), kt(2), mu(2), N);
    
    F.w = [w1; w2];
    F.F = [T1, T2, N];
end