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

%%
%% test if CB modes and matrices are the same with doc.
% read from file
clear
ReadFromCSV;
CB_doc = CB;
Rx_doc = Rx;
xe0_doc = xe0;

% calculate from Mee file
Data;
CriagBamptonReduction;

%% compare CBmodes
eps_Kaa = norm(CB_doc.CBmods.Kaa - CB.CBmods.Kaa) / norm(CB.CBmods.Kaa)
eps_Psi = norm(full(CB_doc.CBmods.Psi - CB.CBmods.Psi)) / norm(full(CB.CBmods.Psi))
eps_Phi = norm(full(CB_doc.CBmods.Phi - CB.CBmods.Phi)) / norm(full(CB.CBmods.Phi))
eps_Max = norm(full(CB_doc.CB_MK.Max - CB.CB_MK.Max)) / norm(full(CB.CB_MK.Max))
eps_Mxx = norm(full(CB_doc.CB_MK.Mxx - CB.CB_MK.Mxx)) / norm(full(CB.CB_MK.Mxx))
eps_Kxx = norm(full(CB_doc.CB_MK.Kxx - CB.CB_MK.Kxx)) / norm(full(CB.CB_MK.Kxx))
eps_Fa = norm(full(CB_doc.CB_F.Fa - CB.CB_F.Fa)) / norm(full(CB.CB_F.Fa))
eps_Fx = norm(full(CB_doc.CB_F.Fx - CB.CB_F.Fx)) / norm(full(CB.CB_F.Fx))

eps_Rx = norm(Rx_doc - Rx) / norm(Rx)
eps_xe0 = norm(xe0_doc - xe0) / norm(xe0)

%% full stuck frq by increasing contact pairs
% clear
% c = jet(5);
% for j = 1:5
%     figure;
%     for i = 4:4:36
%         names = 'full_stuck_modes_CP' + string(i) + '.mat';
%         load(names);
%         plot(i, sqrt(e(j)), 'Color', c(j,:), 'Marker', 'o','MarkerFaceColor', c(j,:)), hold on;
%     end
%     grid on;
%     xlabel('number of contact pairs');
%     ylabel('Omega_' + string(j));
%     title('frequency of mode ' + string(j) + ' VS number of contact pairs');
% end

clear
close all
frqs = zeros(20,9);
for i = 4:4:36
    names = 'full_stuck_modes_CP' + string(i) + '.mat';
    load(names);
    frqs(:, i/4) = e;
end

% figure
% c = jet(5);
% CP = [4:4:36];
% for i = 1:5
%     plot(CP, frqs(i, :), 'Color', c(i,:), 'Marker', 'o','MarkerFaceColor', c(i,:), 'LineStyle', '-'), hold on;
% end
% xlabel('number of contact pairs');
% ylabel('Omega');
% title('frequency of first 5 modes VS number of contact pairs');
% grid on;

% path = pwd;
% cd ../pictures/NewMesh/
c = jet(5);
CP = [4:4:36];

for i = 1:5
    figure
    plot(CP, sqrt(frqs(i, :)), 'Color', c(i,:), 'Marker', 'o','MarkerFaceColor', c(i,:), 'LineStyle', '-'), hold on;

    xlabel('number of contact pairs');
    ylabel('Omega_' + string(i));
    titlename = 'frequency of mode ' + string(i) + ' VS number of contact pairs';
    title(titlename);
    grid on;
    % savefig(titlename);
end

%% plot a1 resonance curve 
for CBmode = 1:1
    figure(1)
    dof_plot = CBmode;
    for i = 4:4:36
        names = 'data/NewMesh/Pe100each_Adof_CP' + string(i) + '_PreDisAll_PreloadFixed.mat';
        load(names);
        plot(Adof(:, 1), Adof(:, dof_plot + 2),'--'), hold on; 
    end
    grid on;
    xlabel('Omega');
    namey = '|dof' + string(dof_plot) + '|';
    ylabel(namey);
    % legend('CP4', 'CP8', 'CP12', 'CP16', 'CP20', 'CP24', 'CP28', 'CP32', 'CP36');
    xlim([4050,4400]);
    % plot([4218,4218], [0, 300], 'k--','DisplayName','Omega4218');
end
%% plot a1 resonance curve - preload act only in first 4 contacts
figure
for i = 4:4:24
    names = 'data/NewMesh/Adof_CP' + string(i) + '.mat';
    load(names);
    plot(Adof(:,1), Adof(:,3)), hold on;
end
grid on;
xlabel('Omega');
ylabel('|a1|');
legend('CP4', 'CP8', 'CP12', 'CP16', 'CP20', 'CP24')

%% the difference of two Rx
Rx = zeros(3*36, 2*8);
j = 1;
for Nx = 8:4:36
    dataname = 'FEM/CP' + string(Nx) + '.mat';
    load(dataname);

    Xcn0 = zeros(3*Nx, 1);
    X00 = [0,0,0.001]';
    for i = 1:Nx
        Xcn0(3 * i - 2:3 * i) = X00;
    end

    FEM.Pe = FEM.Kec * Xcn0 ./ (Nx/4);
    FEM.Pc = FEM.Kcc * Xcn0 ./ (Nx/4);
    xe0 = FEM.Kee \ FEM.Pe; 
    Rx2 = FEM.Kec' * xe0 - FEM.Pc; 


    Rx(1:3*Nx, 2*j - 1) = Rx2;
    j = j + 1;
end

j = 1;
for Nx = 8:4:36
    dataname = 'FEM/CP' + string(Nx) + '.mat';
    load(dataname);

    Xcn0 = zeros(3*Nx, 1);
    X00 = [0,0,0.001]';
    for i = 1:4
        Xcn0(3 * i - 2:3 * i) = X00;
    end

    FEM.Pe = FEM.Kec * Xcn0;
    FEM.Pc = FEM.Kcc * Xcn0;
    xe0 = FEM.Kee \ FEM.Pe; 
    Rx2 = FEM.Kec' * xe0 - FEM.Pc; 

    
    Rx(1:3*Nx, 2*j) = Rx2;
    j = j + 1;
end

%% rename the file
for i = 4:4:36
    oldfilename = 'data/NewMesh/Adof_CP' + string(i) + '.mat';
    newfilename = 'data/NewMesh/Adof_CP' + string(i) + '_PreDisOnly4.mat';
    movefile(oldfilename, newfilename);
end


%% plot property of friction effects

% F vs t
for Nx = 16:4:16

    a = floor(sqrt(Nx));
    while mod(Nx, a) ~= 0
        a = a - 1;
    end
    b = Nx / a;

    paraname = 'data/NewMesh/Data_Omega_4217_Nx_' + string(Nx) + '.mat';
    load(paraname);
    figure
    for i = 1:Nx
        subplot(a,b,i) % Cp1
        plot(para.t, para.Ft(:, 3 * i - 2), 'b-'), hold on;
        plot(para.t, para.Ft(:, 3 * i - 1), 'r-'), hold on;
        muFn = para.fc.mu(1, 1) .* para.Ft(:, 3 * i);
        plot(para.t, muFn, 'k-'), hold on;
        plot(para.t, -muFn, 'k-'), hold on;
    end
end

%%
for i = [1,3,4,5,9]
    disp(i);
end

%% Pc-CP; Xp-CP
% Pe = zeros(36 * 3, 9);
% Pc = zeros(36 * 3, 9);
PeSum = zeros(3, 9); % x, y, z 3-direction in 9 different contacts
PcSum = zeros(3, 9); % t1, t2, n 3-direction in 9 different contacts
Xp = zeros(36 * 3, 9);
figure; % Pc-CP
for i = 4:4:36
    dataname = 'data/NewMesh/Data_Omega_4218_Nx_' + string(i) + '.mat';
    load(dataname);
    % Pe(1:i * 3, i / 4) = para.Pe;
    % Pc(1:i * 3, i / 4) = para.Pc;
    Xp(1:i * 3, i / 4) = para.xp;
    plot(i * ones(i*3, 1), para.Pc, '.'), hold on;

    for k = 1:3
        PeSum(k, i / 4) = sum(abs(para.Pe(k:3:end))) ./ 2;
        PcSum(k, i / 4) = sum(para.Pc(k:3:end));
    end
end
xlabel('CP');
ylabel('Pc');
xticks(0:4:40);
grid on;

% plot Pe sum up in x, y and z directions
figure; 
for i = 1:3
    plot([4:4:36], PeSum(i, :), '-o'), hold on;
end
xlabel('CP');
ylabel('PeSum');
xticks(0:4:40);
legend('x', 'y', 'z');
grid on;

% plot Pc sum up in normal and two tangential directions
figure; 
for i = 1:3
    plot([4:4:36], PcSum(i, :), '-o'), hold on;
end
xlabel('CP');
ylabel('PcSum');
xticks(0:4:40);
legend('T1', 'T2', 'Normal');
grid on;



figure; % xp-CP
for i = 4:4:36
    for j = 3*i-11 : 3*i
        plot([i:4:36], Xp(j, i/4:end), '-o'), hold on;
    end
end
xlabel('CP');
ylabel('xp');
xticks(0:4:40);
grid on;

%% compare the displacements of first 4 CP

linestyles = {'-','--',':','-.'};
markers    = {'o','s','d','^'};

figure; 
for i = 4:4:36
    dataname = 'data/NewMesh/Data_Omega_4218_Nx_' + string(i) + '.mat';
    load(dataname);

    for j = 1:12
        subplot(3,4,j);
        plot(para.t, para.xt(:, 5 + j), 'LineStyle',linestyles{1 + mod(i/4,4)}), hold on;
    end
    
end

for j = 1:12
    subplot(3,4,j);
    grid on;
    legend('CP4', 'CP8', 'CP12', 'CP16', 'CP20', 'CP24', 'CP28', 'CP32', 'CP36');
    xlabel('t');
    ylabel('displacement');
end

%% plot property of friction effects
omega_plot = 4210;

% F vs t
for Nx = 4:4:36

    switch Nx
        case 4
            omega_plot = 4210;
        case 8
            omega_plot = 4214.8;
        case 12
            omega_plot = 4216.2;
        case 16
            omega_plot = 4185;
        case 20 
            omega_plot = 4167.52;
        case 24
            omega_plot = 4175.41;
        case 28
            omega_plot = 4163.58;
        case 32
            omega_plot = 4170;
        case 36
            omega_plot = 4170;
    end
    
    paraname = 'data/NewMesh/Data_PeFixed_Omega_'+ string(omega_plot) + '_Nx_' + string(Nx) + '.mat';
    load(paraname);
    % plot T1, T2 and N forces
    figure

    a = floor(sqrt(Nx));
    while mod(Nx, a) ~= 0
        a = a - 1;
    end
    b = Nx / a;

    for i = 1:Nx
        subplot(a,b,i) % Cp1
        plot(para.t, para.Ft(:, 3 * i - 2), 'b-'), hold on;
        plot(para.t, para.Ft(:, 3 * i - 1), 'r-'), hold on;
        muFn = para.fc.mu(1, 1) .* para.Ft(:, 3 * i);
        plot(para.t, muFn, 'k-'), hold on;
        plot(para.t, -muFn, 'k-'), hold on;
        % legend('T1', 'T2', 'muFn', '-muFn');
        title(i);
    end
end

%% plot x-F
omega_plot = 4210;

% F vs t
for Nx = 36:4:36

    switch Nx
        case 4
            omega_plot = 4210;
        case 8
            omega_plot = 4214.8;
        case 12
            omega_plot = 4216.2;
        case 16
            omega_plot = 4185;
        case 20 
            omega_plot = 4167.52;
        case 24
            omega_plot = 4175.41;
        case 28
            omega_plot = 4163.58;
        case 32
            omega_plot = 4170;
        case 36
            omega_plot = 4170;
    end
    
    paraname = 'data/NewMesh/Pe100eachData_PeFixed_Omega_'+ string(omega_plot) + '_Nx_' + string(Nx) + '.mat';
    load(paraname);
    
    a = floor(sqrt(Nx));
    while mod(Nx, a) ~= 0
        a = a - 1;
    end
    b = Nx / a;
    
    % plot displacements - forces T1
    figure

    for i = 1:Nx
        subplot(a,b,i) % Cp1
        plot(para.xtpxp(:, 3 * i - 2), para.Ft(:, 3 * i - 2), 'b-'), hold on;
        % plot(para.xtpxp(:, 3 * i - 1), para.Ft(:, 3 * i - 1), 'r-'), hold on;
        title(i);
    end

    % plot displacements - forces T2
    figure

    for i = 1:Nx
        subplot(a,b,i) % Cp1
        % plot(para.xtpxp(:, 3 * i - 2), para.Ft(:, 3 * i - 2), 'b-'), hold on;
        plot(para.xtpxp(:, 3 * i - 1), para.Ft(:, 3 * i - 1), 'r-'), hold on;
        title(i);
    end
end

%% plot a1 resonance curve comparison

c = jet(9);
% CP = [4:4:36];
% for i = 1:5
%     plot(CP, frqs(i, :), 'Color', c(i,:), 'Marker', 'o','MarkerFaceColor', c(i,:), 'LineStyle', '-'), hold on;
% end

for CBmode = 1:1
    figure
    dof_plot = CBmode;
    for i = 4:4:36
        names = 'data/NewMesh/Lower_Adof_CP' + string(i) + '_PreDisAll_PreloadFixed.mat';
        load(names);
        plot(Adof(:, 1), Adof(:, dof_plot + 2),'Color', c(i/4,:),'LineStyle','-'), hold on; 
        names = 'data/NewMesh/Pe100each_Adof_CP' + string(i) + '_PreDisAll_PreloadFixed.mat';
        load(names);
        plot(Adof(:, 1), Adof(:, dof_plot + 2),'Color', c(i/4,:),'LineStyle','--'), hold on; 
    end
    grid on;
    xlabel('Omega');
    namey = '|dof' + string(dof_plot) + '|';
    ylabel(namey);
    legend('CP4','CP4', 'CP8', 'CP8', 'CP12','CP12', 'CP16','CP16', 'CP20','CP20', 'CP24','CP24', 'CP28','CP28', 'CP32','CP32', 'CP36', 'CP36');
    xlim([4050,4400]);
    % plot([4218,4218], [0, 300], 'k--','DisplayName','Omega4218');
end