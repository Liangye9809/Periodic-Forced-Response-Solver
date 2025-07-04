function F = gg(x1, x2, x3, pfunc)
    fc = pfunc.fc;
    
    kn = fc.kn; % stiffness in normal direction
    xn0 = fc.xn0; % initial gap between this contact point
    mu = fc.mu; % friciton coefficient of coulomb friction
    kt = fc.kt; % stiffness of 2 tangential direction in column vector
    w = fc.w;

    