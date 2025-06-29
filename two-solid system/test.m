t = [0:100] * (2*pi/100);
y1 = exp(t/2);
plot(t,y1);
grid on;

%% test handle variables
clear testfunc
testfunc.fc = CoulombFrictionParas(Coulombstruct);
testfunc.fc.w = CoulombFrictionW([1;1], Nx);




function F = testgf(x, p)
    fc = p.fc;

    kn = fc.kn; % 
    xn0 = fc.xn0; 
    mu = fc.mu; % 
    kt = fc.kt; % 
    w = fc.w;

    w = w + 1;
    
    fc.w = w;

    F = w;
end
%%
tic
testF = g(rand(256,17), 0);
toc
%%
clc
disp(['dt   ','nsteps   ','max(dtA) ','max(expm(dtA))    ', 'max(abs(eig(expm(dtA))))    ', 'max(abs(eigEuler))  ', 'max(abs(eigRK2))', 'ImpEular']);
clear J
J = finite_diff_jac(@(x) g(x, params.func).F, xp);
J = [zeros(5,17);
      zeros(12,5), J];
J = -M \ J;
Jg = [zeros(17,17), zeros(17,17);
      J, zeros(17,17)];
for i  = 0:8
    nsteps = 100 * 10^i;
    % nsteps = 100 * 80;
    testdt = T / nsteps;
    testdtA = testdt * A;
    R = expm(testdtA) + testdt * expm(testdtA) * Jg;
    e = eig(R);
    eA = eig((eye(34) - testdt * A));
    RK2R = testdt * A + testdt^2 / 2 * A^2 + eye(34);
    ImpR = inv(eye(34) - testdt*(A + Jg));
    eImpR = eig(ImpR);
    maxe = max(abs(e));
    maxeE = max(abs(eA));
    maxeRK2 = max(abs(eig(RK2R)));
    disp([testdt, nsteps, max(max(abs(testdtA))), max(max(abs(R))), maxe, maxeE, maxeRK2, max(abs(eImpR))]);
end
%% 
t = -1:0.001:1;
y = sqrt(abs(t));
plot(t,y)
%%
a_1 = 11/6;
a0 = -3;
a1 = 3/2;
a2 = -1/3;
f1 = a_1 + a0 + a1 + a2;
f2 = a0 + 2*a1 + 3*a2;
f3 = a0 + 4*a1 + 9*a2;
f4 = a0 + 8*a1 + 27*a2;
%%
for i  = 8:8
    nsteps = 14500;
    testdt = T / nsteps;
    R = eye(34) + (testdt/2)*A;
    format short g
    disp([testdt, nsteps, max(max(abs(R))), max(max(abs(eig(R))))]);
end