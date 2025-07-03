clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% x¨ + x' + x + 10x^3 = cos(omega*t) %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = Phi_ty(t, y, omega)

    G = [0; 10*y(2)^3];
    F = [0; cos(omega * t)];
    phi = - G + F;

end

A = [0, 1;
     -1, -1];

omega_0 = 1;
omega_end = 2;
nstep = 30;
nrelax = 10;
nomega = 100;
y = zeros(2, 1+nstep*nrelax, nomega+1);
ye = zeros(2, 1+nstep*nrelax, nomega+1);
t(1,1) = 0;
i = 1;
omega_cont = [];
domega = (omega_end - omega_0) / nomega;

tic;
for omega = omega_0:domega:omega_end
    omega_cont = [omega_cont, omega];
    T = 2*pi / omega;
    dt = T / nstep;
    R = expm(dt*A); 
    Re = (eye(2) + dt*A);
    for k = 1:nrelax*nstep
        y(:, k+1, i) = R * (y(:, k, i) + dt * Phi_ty(t(i, k), y(:, k, i), omega));
        ye(:, k+1, i) = Re * ye(:, k, i) + dt * Phi_ty(t(i, k), y(:, k, i), omega);
        t(i, k+1) = t(i, k) + dt;
    end
    disp('omega ' + string(omega) + ' a1(end) ' + string(y(1, end, i)));
    if omega >= omega_end
        break;
    end
    i = i + 1;
    t(i, 1) = 0;
    y(:, 1, i) = y(:, end, i-1);
    ye(:, 1, i) = ye(:, end, i-1);
end

toc;

figure
plot(y(1,:,1),'b'), hold on;
plot(ye(1,:,1), 'r'), hold on
%%
clear Amax
for j = 1:size(y,3) % omega number
    Amax(:, j) = max(abs(y(:,end-10*nstep:end,j)), [], 2);
end
% omega_cont = omega_0:domega:omega_end;
figure
plot(omega_cont, Amax(1,:),'ko');
grid on;
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
grid on;