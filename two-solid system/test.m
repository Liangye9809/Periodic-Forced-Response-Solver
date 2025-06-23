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
disp(['dt   ','nsteps   ','max(dtA) ','max(exp(dtA))    ', 'max(abs(eig(exp(dtA))))    ', 'max(abs(eigEuler))  ', 'max(abs(eigRK2))']);
for i  = 0:10
    nsteps = 100 * 10^i;
    testdt = T / nsteps;
    testdtA = testdt * A;
    testexpdtA = exp(testdtA);
    e = eig(testexpdtA);
    eA = eig((eye(34) - testdt * A));
    RK2R = testdt * A + testdt^2 / 2 * A^2 + eye(34);
    maxe = max(abs(e));
    maxeE = max(abs(eA));
    maxeRK2 = max(abs(eig(RK2R)));
    disp([testdt, nsteps, max(max(abs(testdtA))), max(max(abs(testexpdtA))), maxe, maxeE, maxeRK2]);
end
%% 
t = -1:0.001:1;
y = sqrt(abs(t));
plot(t,y)