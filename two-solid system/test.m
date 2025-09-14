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
% clc
disp(['dt   ','nsteps   ','']);
clear J
J = finite_diff_jac(@(x) gf(x, params.func).F, xp);
J = [zeros(5,17);
      zeros(12,5), J];
J = -M \ J;
Jg = [zeros(17,17), zeros(17,17);
      J, zeros(17,17)];
T = 2*pi / omega_0;
for i  = 0:8
    nsteps = 100 * 10^i;
    dt = T / nsteps;

    % explicit Euler exponential numerical method
    dtA = dt * A;
    R_ExpEulxpm = expm(dtA) + dt * expm(dtA) * Jg;
    e_ExpEulxpm = eig(R_ExpEulxpm);
    maxe_ExpEulxpm = max(abs(e_ExpEulxpm));

    % % RK2 exponential numerical method
    % RK2R = dt * A + dt^2 / 2 * A^2 + eye(34);
    % maxeRK2 = max(abs(eig(RK2R)));


    % Implicit Euler Method
    R_ImpEulxpm = inv(eye(34) - dt*(A + Jg));
    e_ImpEulxpm = eig(R_ImpEulxpm);
    maxe_ImpEulxpm = max(abs(e_ImpEulxpm));
    
    
    disp([dt, nsteps, maxe_ExpEulxpm, max(abs(e_ImpEulxpm))]);
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
    dt = T / nsteps;
    R_ExpEulxpm = eye(34) + (dt/2)*A;
    format short g
    disp([dt, nsteps, max(max(abs(R_ExpEulxpm))), max(max(abs(eig(R_ExpEulxpm))))]);
end

%% test g only
clear
load("ptest.mat");
N = 2^3;
Nx = 4;
xtest = zeros(N, 3 * Nx);
fc.kn = ptest.fc.kn;
fc.xn0 = ptest.fc.xn0;
fc.mu = ptest.fc.mu;
fc.kt = ptest.fc.kt;
fc.w = ptest.fc.w;
for i = 1:N
    for j = 1:(3 * Nx)
        xtest(i, j) = 12 * (i -1) + j;
    end
end
profile on
tic;
for i = 1:1
    F = g_mex(xtest, fc);
end
toc;
profile off
p = profile('info');
profview(0,p);

%%
clear
load("fc_row.mat");
N = 2^8;
Nx = 4;
xtest = rand(N, 3 * Nx);

tic;
for i = 1:1000
    % [F, w] = g_mex(xtest, fc);
    F = Getgx(xtest, fc);
    % F = SeletF(g_mex(xtest, fc));
end
toc;
