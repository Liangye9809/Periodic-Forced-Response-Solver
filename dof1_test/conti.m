clear
clc
close all
%% Implicit Euler Method
tic;
K = 1;
M = 1;
C = 1;
F0 = 1;
A = [  0,  1;
     -M\K, -M\C];
F = M \ (2 * F0);
omega_0 = 0.5;
omega_end = 1.2;
y0 = zeros(2, 1);
y(:,1,1) = y0;
t(1,1) = 0;
i = 1;
omega_cont = [];
domega = (omega_end - omega_0) / 50;
for omega = omega_0:domega:omega_end
    omega_cont = [omega_cont, omega];
    T = 2*pi / omega;
    nstep = 300;
    nrelax = 100;
    dt = T / nstep;
    R = eye(2) - dt*A; % Implicit Euler Method
    for k = 1:nrelax*nstep
        y(:, k+1, i) = R \ (y(:, k, i) + dt * (F * sin(omega * (dt + t(i, k)))));
        t(i, k+1) = t(i, k) + dt;
    end
    if omega >= omega_end
        break;
    end
    i = i + 1;
    t(i, 1) = 0;
    y(:, 1, i) = y(:, end, i-1);
end

toc;
%%
clear Amax
for j = 1:size(y,3) % omega number
    Amax(:, j) = max(abs(y(:,end-10*nstep:end,j)), [], 2);
end
% omega_cont = omega_0:domega:omega_end;
figure
plot(omega_cont, Amax(1,:),'ko');
%%
figure;
i = 1;
tt = 0;
for omega = omega_0:domega:omega_end
    T = 2 * pi / omega;
    dt = T / nstep;
    tt = tt(end) + dt*[0:nstep*nrelax];
    plot(tt, y(1,:,i)), hold on;
    i = i + 1;
end