% function F = gf(x, p) % x is a 3dof column vector, represent 2 tangential directions and 1 normal direction
%     kn = p.kn; % stiffness in normal direction
%     xn0 = p.xn0; % initial gap between this contact point
%     mu = p.mu; % friciton coefficient of coulomb friction
%     kt = p.kt; % stiffness of 2 tangential direction in column vector
%     % w = p. w; % slip displacement of 2 tangential direction in column vector
%     N = NormalForces(x(3), kn, xn0);
%     [T1, p.w(1)] = TangentialForce(x(1), p.w(1), kt(1), mu(1), N);
%     [T2, p.w(2)] = TangentialForce(x(2), p.w(2), kt(2), mu(2), N);
%     F = [T1; T2; N];
% end

function F = gf(x, p) % x is a 3*Np dof row vector, represent 2*Np tangential directions and 1*Np normal direction

    kn = p.fc.kn; % stiffness in normal direction
    xn0 = p.fc.xn0; % initial gap between this contact point
    mu = p.fc.mu; % friciton coefficient of coulomb friction
    kt = p.fc.kt; % stiffness of 2 tangential direction in column vector
    Np = size(kn, 2);
    FN = zeros(1, Np);
    T1 = zeros(1, Np);
    T2 = zeros(1, Np);
    for i = 1:Np
        FN(i) = NormalForces(x(3*(i-1)+3), kn(i), xn0(i));
        [T1(i), p.fc.w(1, i)] = TangentialForces(x(3*(i-1)+1), p.fc.w(1, i), kt(1, i), mu(1, i), FN(i));
        [T2(i), p.fc.w(2, i)] = TangentialForces(x(3*(i-1)+2), p.fc.w(2, i), kt(2, i), mu(2, i), FN(i));
    end
    F = [T1; T2; FN];
    F = F(:);
end
