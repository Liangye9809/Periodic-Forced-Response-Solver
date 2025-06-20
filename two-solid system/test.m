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
disp(['dt   ','nsteps   ','max(dtA) ','max(exp(dtA))']);
for i  = 0:5
    nsteps = 100 * 10^i;
    testdt = T / nsteps;
    testdtA = testdt * A;
    testexpdtA = exp(testdtA);
    disp([testdt, nsteps, max(max(abs(testdtA))), max(max(abs(testexpdtA)))]);
end
%% 
t = -1:0.001:1;
y = sqrt(abs(t));
plot(t,y)